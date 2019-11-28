# animate_interactive.py
#
# Assemble a 2D animation from a 3D array.
#
# Author: S. Candelaresi (iomsn1@gmail.com).
"""
Versatile interactive plotting routine for slices and tomography.
"""


def animate_interactive(data, t=None, dim_order=(0, 1, 2),
                        fps=10.0, title=None, xlabel='x', ylabel='y',
                        font_size=24, color_bar=0, colorbar_label=None,
                        sloppy=True, fancy=False,
                        range_min=None, range_max=None, extent=[-1, 1, -1, 1],
                        shade=False, azdeg=0, altdeg=65,
                        arrowsX=None, arrowsY=None, arrows_resX=10, arrows_resY=10,
                        arrows_pivot='mid', arrows_width=0.002, arrows_scale=5,
                        arrows_color='black', plot_arrows_grid=False,
                        movie_file=None, bitrate=1800, keep_images=False,
                        figsize=(8, 7), dpi=300,
                        **kwimshow):
    """
    Assemble a 2D animation from a 3D array.

    call signature::

    animate_interactive(data, t=None, dim_order=(0, 1, 2),
                        fps=10.0, title=None, xlabel='x', ylabel='y',
                        font_size=24, color_bar=0, colorbar_label=None,
                        sloppy=True, fancy=False,
                        range_min=None, range_max=None, extent=[-1, 1, -1, 1],
                        shade=False, azdeg=0, altdeg=65,
                        arrowsX=None, arrowsY=None, arrows_resX=10, arrows_resY=10,
                        arrows_pivot='mid', arrows_width=0.002, arrows_scale=5,
                        arrows_color='black', plot_arrows_grid=False,
                        movie_file=None, bitrate=1800, keep_images=False,
                        figsize=(8, 7), dpi=300,
                        **kwimshow)

    Assemble a 2D animation from a 3D array. *data* has to be a 3D array of
    shape [nt, nx, ny] and who's time index has the same dimension as *t*.
    The time index of *data* as well as its x and y indices can be changed
    via *dim_order*.

    Keyword arguments:

    *dim_order*:
      Ordering of the dimensions in the data array (t, x, y).

    *fps*:
      Frames per second of the animation.

    *title*:
      Title of the plot.

    *xlabel*:
      Label of the x-axis.

    *ylabel*:
      Label of the y-axis.

    *font_size*:
      Font size of the title, x and y label.
      The size of the x- and y-ticks is 0.5*font_size and the colorbar ticks'
      font size is 0.5*font_size.

    *color_bar*: [ 0 | 1 ]
      Determines how the colorbar changes:
      (0 - no cahnge; 1 - adapt extreme values).

    *colorbar_label*:
      Label of the color bar.

    *sloppy*: [ True | False ]
      If True the update of the plot lags one frame behind. This speeds up the
      plotting.

    *fancy*: [ True | False ]
      Use fancy font style.

    *range_min*, *range_max*:
      Range of the colortable.

    *extent*: [ None | (left, right, bottom, top) ]
      Limits for the axes (domain).

    *shade*: [ False | True ]
      If True plot a shaded relief instead of the usual colormap.
      Note that with this option cmap has to be specified like
      cmap = plt.cm.hot instead of cmap = 'hot'. Shading cannot
      be used with the color_bar = 0 option.

    *azdeg*, *altdeg*:
      Azimuth and altitude of the light source for the shading.

    *arrowsX*:
      Data containing the x-component of the arrows.

    *arrowsY*:
      Data containing the y-component of the arrows.

    *arrows_resXY*:
      Plot every arrows_resXY arrow in x and y.

    *arrows_pivot*: [ 'tail' | 'middle' | 'tip' ]
      The part of the arrow that is used as pivot point.

    *arrows_width*:
      Width of the arrows.

    *arrows_scale*:
      Scaling of the arrows.

    *arrows_color*:
      Color of the arrows.

    *plot_arrows_grid*: [ False | True ]
      If 'True' the grid where the arrows are aligned to is shown.

    *movie_file*: [ None | string ]
      The movie file where the animation should be saved to.
      If 'None' no movie file is written. Requires 'ffmpeg' to be installed.

    *bitrate*:
      Bitrate of the movie file. Set to higher value for higher quality.

    *keep_images*: [ False | True ]
      If 'True' the images for the movie creation are not deleted.

    *figsize*:
      Size of the figure in inches.

    *dpi*:
      Dots per inch of the frame.

    **kwimshow:
      Remaining arguments are identical to those of pylab.imshow. Refer to that help.
    """

    try:
        import thread
    except:
        import _thread as thread

    # We need to define these variables as globals, as they are being used
    # by various functions.

    global time_step, time_slider, pause
    global fig, axes, image, colorbar, arrows, manager, n_times, movie_files
    global rgb, plot_arrows

    if title is None:
        title = ''


    def plot_frame():
        """
        Plot the current frame.
        """

        global time_step, axes, colorbar, arrows, manager, rgb

        # Define the plot title.
        if not movie_file is None:
            axes.set_title(title + r'$\quad$' + r'$t={0:.4e}$'.format(t[time_step]),
                           fontsize=font_size)

        # Update the image data.
        if not shade:
            image.set_data(data[time_step, :, :])
        else:
            image.set_data(rgb[time_step, :, :, :])

        # Update the colorbar.
        if color_bar == 0:
            pass
        if color_bar == 1:
            colorbar.set_clim(vmin=data[time_step, :, :].min(),
                              vmax=data[time_step, :, :].max())
            colorbar.draw_all()

        # Update the arrows data.
        if plot_arrows:
            arrows.set_UVC(U=arrowsX[time_step, ::arrows_resX, ::arrows_resY],
                           V=arrowsY[time_step, ::arrows_resX, ::arrows_resY])

        if not sloppy or (not movie_file is None):
            manager.canvas.draw()


    def play(thread_name):
        """
        Play the movie.
        """

        import time
        global time_step, time_slider, pause, fig, axes, n_times, movie_files

        pause = False
        while (time_step < n_times) and (not pause):
            # Write the image files for the movie.
            if not movie_file is None:
                plot_frame()
                frame_name = '{0}{1:06}.png'.format(movie_file, time_step)
                fig.savefig(frame_name, dpi=dpi)
                movie_files.append(frame_name)
            else:
                time_start = time.clock()
                time_slider.set_val(t[time_step])
                # Wait for the next frame (fps).
                while (time.clock() - time_start < 1.0/fps):
                    pass
            time_step += 1
        time_step -= 1


    def play_thread(event):
        """
        Call the play function as a separate thread (for GUI).
        """

        global pause

        if pause:
            try:
                thread.start_new_thread(play, ("play_thread", ))
            except:
                print("Error: unable to start play thread.")


    def pausing(event):
        global pause

        pause = True


    def reverse(event):
        global time_step, time_slider

        time_step -= 1
        if time_step < 0:
            time_step = 0
        # Plot the frame and update the time slider.
        time_slider.set_val(t[time_step])


    def forward(event):
        global time_step, time_slider

        time_step += 1
        if time_step > len(t)-1:
            time_step = len(t)-1
        # Plot the frame and update the time slider.
        time_slider.set_val(t[time_step])


    import numpy as np
    import pylab as plt

    pause = True
    plot_arrows = False

    # Check if the data has the right dimensions.
    if (data.ndim != 3 and data.ndim != 4):
        print("Error: data dimensions are invalid: {0} instead of 3.".format(data.ndim))
        return -1

    # Transpose the data according to dim_order.
    unordered_data = data
    data = np.transpose(unordered_data, dim_order)
    del(unordered_data)

    # Check if arrows should be plotted.
    if not(arrowsX is None) and not(arrowsY is None):
        if (isinstance(arrowsX, np.ndarray) and isinstance(arrowsY, np.ndarray)):
            if arrowsX.ndim == 3:
                # Transpose the data according to dim_order.
                unordered_data = arrowsX
                arrowsX = np.transpose(unordered_data, dim_order)
                del(unordered_data)
            if arrowsY.ndim == 3:
                # Transpose the data according to dim_order.
                unordered_data = arrowsY
                arrowsY = np.transpose(unordered_data, dim_order)
                unordered_data = []

                # Check if the dimensions of the arrow arrays match each other.
                if arrowsX.shape != arrowsY.shape:
                    print("Error: dimensions of arrowX do not match with dimensions of arrowY.")
                    return -1
                else:
                    plot_arrows = True
        else:
            print("Warning: arrowsX and/or arrowsY are of invalid type.")

    # Check if time array has the right length.
    n_times = len(t)
    if n_times != data.shape[0]:
        print("Error: length of time array does not match length of data array.")
        return -1
    if plot_arrows:
        if (n_times != arrowsX.shape[0]) or (n_times != arrowsY.shape[0]):
            print("error: length of time array does not match length of arrows array.")
            return -1

    # Check if fps is positive.
    if fps < 0:
        print("Error: fps is not positive, fps = {0}.".format(fps))
        return -1

    # Determine the size of the data array.
    nX = data.shape[1]
    nY = data.shape[2]

    # Determine the minimum and maximum values of the data set.
    if not range_min:
        range_min = np.min(data)
    if not range_max:
        range_max = np.max(data)

    # Setup the plot.
    if fancy:
        plt.rc('text', usetex=True)
        plt.rc('font', family='arial')
    else:
        plt.rc('text', usetex=False)
        plt.rc('font', family='sans')
    if not movie_file is None:
        fig = plt.figure(figsize=figsize)
        axes = plt.axes([0.15, 0.1, .70, .85])
    else:
        fig = plt.figure(figsize=figsize)
        axes = plt.axes([0.1, 0.3, .80, .65])

    # Set up canvas of the plot.
    axes.set_title(title, fontsize=font_size)
    axes.set_xlabel(xlabel, fontsize=font_size)
    axes.set_ylabel(ylabel, fontsize=font_size)
    plt.xticks(fontsize=0.5*font_size)
    plt.yticks(fontsize=0.5*font_size)
    if shade:
        plane = np.zeros([nX, nY, 3])
    else:
        plane = np.zeros([nX, nY])

    # Apply shading.
    if shade:
        from matplotlib.colors import LightSource

        ls = LightSource(azdeg=azdeg, altdeg=altdeg)
        rgb = []
        # Shading can be only used with color_bar=1 or color_bar=2 at the moment.
        if color_bar == 0:
            color_bar = 1
        # Check if colormap is set, if not set it to 'copper'.
        if 'cmap' not in kwimshow.keys():
            kwimshow['cmap'] = plt.cm.copper
        for i in range(data.shape[0]):
            tmp = ls.shade(data[i, :, :], kwimshow['cmap'])
            rgb.append(tmp.tolist())
        rgb = np.array(rgb)
        del(tmp)

    # Calibrate the displayed colors for the data range.
    image = axes.imshow(plane, vmin=range_min, vmax=range_max, origin='lower',
                        extent=extent, **kwimshow)
    colorbar = fig.colorbar(image)
    colorbar.set_label(colorbar_label, fontsize=font_size, labelpad=10)

    # Change the font size of the colorbar's ytickslabels.
    cbytick_obj = plt.getp(colorbar.ax.axes, 'yticklabels')
    plt.setp(cbytick_obj, fontsize=0.5*font_size)

    # Plot the arrows.
    if plot_arrows:
        # Prepare the mesh grid where the arrows will be drawn.
        arrow_grid = np.meshgrid(np.arange(extent[0], extent[1],
                                           float(extent[1]-extent[0])*arrows_resX/(data.shape[2]-1)),
                                 np.arange(extent[2], extent[3],
                                           float(extent[3]-extent[2])*arrows_resY/(data.shape[1]-1)))
        arrows = axes.quiver(arrow_grid[0], arrow_grid[1],
                             arrowsX[0, ::arrows_resX, ::arrows_resY],
                             arrowsY[0, ::arrows_resX, ::arrows_resY],
                             units='width', pivot=arrows_pivot, width=arrows_width,
                             scale=arrows_scale, color=arrows_color)
        # Plot the grid for the arrows.
        if plot_arrows_grid:
            axes.plot(arrow_grid[0], arrow_grid[1], 'k.')

    # For real-time image display.
    if (not sloppy) or (not movie_file is None):
        manager = plt.get_current_fig_manager()
        manager.show()

    time_step = 0
    if not movie_file is None:
        import os

        movie_files = []

        # Start the animation.
        play('no_thread')

        # Write the movie file.
        ffmpeg_command = "ffmpeg -r {0} -i {1}%6d.png -vcodec mpeg4 -b:v {2} -q:v 0 {3}.avi".format(fps, movie_file, bitrate, movie_file)
        os.system(ffmpeg_command)
        # Clean up the image files.
        if not keep_images:
            print("Cleaning up files.")
            for fname in movie_files:
                os.remove(fname)
    else:
        # Set up the gui.
        plt.ion()
        plt.subplots_adjust(bottom=0.2)

#        axes_play = plt.axes([0.1, 0.05, 0.15, 0.05])
#        button_play = plt.Button(axes_play, 'play', color='lightgoldenrodyellow',
#                                 hovercolor='0.975')
#        button_play.on_clicked(play_thread)

#        axes_pause = plt.axes([0.3, 0.05, 0.15, 0.05])
#        button_pause = plt.Button(axes_pause, 'pause', color='lightgoldenrodyellow',
#                                  hovercolor='0.975')
#        button_pause.on_clicked(pausing)

#        axes_reverse = plt.axes([0.5, 0.05, 0.15, 0.05])
        axes_reverse = plt.axes([0.1, 0.05, 0.3, 0.05])
        button_reverse = plt.Button(axes_reverse, 'reverse', color='lightgoldenrodyellow',
                                    hovercolor='0.975')
        button_reverse.on_clicked(reverse)

#        axes_forward = plt.axes([0.7, 0.05, 0.15, 0.05])
        axes_forward = plt.axes([0.5, 0.05, 0.3, 0.05])
        button_forward = plt.Button(axes_forward, 'forward', color='lightgoldenrodyellow',
                                    hovercolor='0.975')
        button_forward.on_clicked(forward)

        # Create the time slider.
        time_slider_axes = plt.axes([0.2, 0.12, 0.6, 0.03], facecolor='lightgoldenrodyellow')
        time_slider = plt.Slider(time_slider_axes, 'time', t[0], t[-1], valinit=t[0])
        def update(val):
            global time_step
            # Find the closest time step to the slider time value.
            for i in range(len(t)):
                if t[i] < time_slider.val:
                    time_step = i
            if (time_step != len(t)-1):
                if (t[time_step+1] - time_slider.val) < (time_slider.val - t[time_step]):
                    time_step += 1
            plot_frame()
        time_slider.on_changed(update)

        plt.show()

    print("done")

#    return button_play, button_pause, button_reverse, button_forward
    return button_reverse, button_forward
