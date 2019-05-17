# rvid_box.py
#
# Generate 3D visualisation from 4x 2D video slices.
#
# Author: S. Long (shiting.long@aalto.fi) and F. Gent (fred.gent.ncl@gmail.com)
"""
This script contains three functions, whereas plot_box is the only function
that needs to be called by the user.
One can call plot_box from a Python interpreter after importing rvid_box.py.
An example call: plot_box(field='lnTT', slice=50, datadir=<path>, unit='temperature').
The output image will be saved to directory visu/images/.

The details of the three functions are listed as follows:

1) plot_box : the main function to plot cube objects by assembling slice_xy, slice_xy2, slice_xz and slice_yz together,
it calls read_slices first then plot.

Arguments:

field --- a string variable, indicating which variable to slice.

datadir --- a string variable, indicating the path to the data directory.

proc --- an integer giving the processor to read a slice from.

scale --- a string variable, the scale of the colorbar; default is regular scale,
    can be changed to logarithmic by setting it 'log' (the base number is 10).

quiet --- a boolean variable, print debug information.

slice --- the sequential integer number of a slice in given period (starting from 0), default is the first slice (slice 0),
    can be set to '-1' to generate all slices in given period.

colorscale --- a string variable,  color scale for the plot. Default is 'Hot'.
    Plotly supports most matplotlib colormaps and you can customize a colorscale on your own.
    See https://plot.ly/python/colorscales/ for more information.

unit --- a string variable, define the unit of the colorbar, shown on its title.

axestitle --- a three-tuple indicating the title of x, y and z axis. It is 'x', 'y', 'z' by default.

viewpoint --- a three-tuple indicating the viewpoint of the center of the box (default is (0, 0, 0),
    you can change the center in the layout of the figure in function plot()).

imageformat --- a string variable, indicating the format of the output image. Set to 'png' by default.
    It supports png, jpeg, webp, svg and pdf format.

2) pcn.read.slices: read the four slices
3) plot: plot cube objects

"""


import pencilnew as pcn
import numpy as np
import os
import plotly.graph_objs as go
import plotly.io as pio

def plot(
         #field (or multiple todo) and 4 surfaces (extra todo) 
         slice_obj, fields, xyzplane,
         itt, it, quiet=True,  
         #yz[i], xy[i], xz[i], xy2[i], 
         #save output
         figdir='./images/', imageformat='.png',
         #set color parameters
         norm='linear', colorscale='RdBu',
         color_levels=None,
         color_range=None,
         #locate the box and axes
         viewpoint=(-1.35, -2.1, 0.5), offset=5., margin=(20,20,30,0),
         autosize=False, image_dim=(800,500),
         #handle axes properties
         visxyz=[True,True,True], axestitle=('x','y','z'), xyz=None,
         #color bar properties
         cbar_label='', cbar_loc=1., cbar_label_pos='right', 
         #add text for time stamp
         time=0, textxy=(0,0), unit='', isd=3, fontsize=25,
        ):

    """
    Arguments:
        yzslices --- a numpy array, data of a yz slice.
        xyslices --- a numpy array, data of a xy slice.
        xzslices --- a numpy array, data of a xz slice.
        xy2slices --- a numpy array, data of a xy2 slice.
        scale --- a string variable, scale of the colorbar; default is regular scale,
            can be changed to logarithmic by setting it 'log' (the base number is 10).
        slice --- the sequential integer number of a slice in given period (starting from 0),
            default is the first slice (slice 0), can be set to '-1' to generate all slices in given period.
        colorscale --- a string variable,  color scale for the plot. Default is 'Hot'.
            Plotly supports most Matplotlib colormaps and you can customize a colorscale on your own.
            See https://plot.ly/python/colorscales/ for more information.
        unit --- a string variable, define the unit of the colorbar, shown on its title.
        axestitle --- a three-tuple indicating the title of x, y and z axis. It is 'x', 'y', 'z' by default.
        viewpoint --- a three-tuple indicating the viewpoint of the center of the box (default is (0, 0, 0),
            you can change the center in the layout of the figure in function plot()).
        field --- a string variable, which variable to slice.
        time --- a double variable, the simulation time when the data to be plotted is generated.
        xyz --- a three-tuple indicating the value range of the three axes.
        offset --- a double variable, the offset of xy2 slice below the box.
            It is set to 2 by default, meaning the xy2 slice is located below the box by half of the height of the box.
        imageformat --- a string variable, indicating the format of the output image. Set to 'png' by default.
            It supports png, jpeg, webp, svg and pdf format.
    """
    """
    options = Options()
    options.add_argument("--headless")
    # options.add_argument( "--screenshot test.jpg http://google.com/" )
    driver = webdriver.Firefox(options=options, executable_path="/homeappl/home/longs1/tmp/geckodriver")
    driver.set_window_size(1000, 500)
    time.sleep(5)
    driver.get('file:///homeappl/home/longs1/appl_taito/pencil-code/python/pencilnew/visu/temp-plot.html')
    driver.save_screenshot('test.png')
    #driver.close()
    #imgkit.from_file('test.html', 'out.png')
    #pio.write_image(fig, "test", format="png")
    """


    if not quiet:
        print('Printing t={:.2g} box plot'.format(time))
    for field in fields:
        cmin,cmax=1e38,-1e38
        #set color limits based on time series or single snapshot
        if color_levels=='common':
            for key in xyzplane:
                globals()[key+'slice']=slice_obj.__getattribute__(key).__getattribute__(field)
                cmax = max(cmax,globals()[key+'slice'][it].max())
                cmin = min(cmin,globals()[key+'slice'][it].min())
        else:
            for key in xyzplane:
                globals()[key+'slice']=slice_obj.__getattribute__(key).__getattribute__(field)
                cmax = max(cmax,globals()[key+'slice'][itt].max())
                cmin = min(cmin,globals()[key+'slice'][itt].min())

        height= globals()['xzslice'][itt].shape[0]
        width = globals()['xyslice'][itt].shape[1]
        
        if xyz is None:
            xx = np.linspace(0, width - 1, width)
            yy = np.linspace(0, width - 1, width)
            zz = np.linspace(0, height - 1, height)
            x_z, z_z = np.meshgrid(xx, zz)
            x_y, y_y = np.meshgrid(xx, yy)
        else:
            x_z, z_z = np.meshgrid(xyz[0], xyz[2])
            x_y, y_y = np.meshgrid(xyz[0], xyz[1])

        if norm == 'log':
            # field argu
            if 'ln' in field:
                z1 = globals()['xyslice' ][itt]
                y1 = globals()['xzslice' ][itt]
                x1 = globals()['yzslice' ][itt]
                z2 = globals()['xy2slice'][itt]
            else:
                z1 = np.log10(globals()['xyslice' ][itt])
                y1 = np.log10(globals()['xzslice' ][itt])
                x1 = np.log10(globals()['yzslice' ][itt])
                z2 = np.log10(globals()['xy2slice'][itt])
                cmax=np.log10(cmax)
                cmin=np.log10(cmin)
        elif norm == 'linear': 
            z1 = globals()['xyslice' ][itt]
            y1 = globals()['xzslice' ][itt]
            x1 = globals()['yzslice' ][itt]
            z2 = globals()['xy2slice'][itt]
            cmax = max(-cmin, cmax)
            cmin = -cmax
        else: 
            print("WARNING: 'norm' undefined, applying 'linear'")
            z1 = globals()['xyslice' ][itt]
            y1 = globals()['xzslice' ][itt]
            x1 = globals()['yzslice' ][itt]
            z2 = globals()['xy2slice'][itt]
            cmax = max(-cmin, cmax)
            cmin = -cmax

        # set offsets of four slices for placing them in correct positions
        if xyz is None:
            z1_offset = -height / 2  * np.ones(z1.shape)
            y1_offset =           0  * np.ones(y1.shape)
            x1_offset =           0  * np.ones(x1.shape)
            z2_offset = (height - 1) * np.ones(z2.shape)
        else:
            z1_offset =(xyz[2].min() - (xyz[2].max() - xyz[2].min())/
                            offset)  * np.ones(z1.shape)
            y1_offset = xyz[1].min() * np.ones(y1.shape)
            x1_offset = xyz[0].min() * np.ones(x1.shape)
            z2_offset = xyz[2].max() * np.ones(z2.shape)

        # projection in the z-direction
        proj_z = lambda x, y, z: z

        # for projection slice xy
        colorsurfz1 = proj_z(x_z, z_z, z1.tolist())
        colorsurfz2 = proj_z(x_z, z_z, z2.tolist())

        # for projection slice xz
        colorsurfy1 = proj_z(x_y, y_y, y1.tolist())

        # for projection slice yz
        colorsurfx1 = proj_z(x_y, y_y, x1.tolist())

        # plot slices to surfaces
        trace_y1 = go.Surface(z=list(z_z),
                              x=list(x_z),
                              y=list(y1_offset),
                              showscale=True,
                              colorscale=colorscale,
                              surfacecolor=colorsurfy1,
                              cmin=cmin,
                              cmax=cmax,
                              colorbar=dict(
                                  title=cbar_label,
                                  x=cbar_loc,
                              )
                              )
        trace_y1.colorbar.title.side=cbar_label_pos
        trace_y1.colorbar.title.font.size=fontsize
 
        trace_x1 = go.Surface(z=list(z_z),
                              x=list(x1_offset),
                              y=list(x_z),
                              showscale=False,
                              surfacecolor=colorsurfx1,
                              colorscale=colorscale,
                              cmin=cmin,
                              cmax=cmax,
                              colorbar=dict(
                              #    title=cbar_label,
                                  x=cbar_loc,
                              )
                              )
        trace_z1 = go.Surface(z=list(z1_offset),
                              x=list(x_y),
                              y=list(y_y),
                              showscale=False,
                              surfacecolor=colorsurfz1,
                              colorscale=colorscale,
                              cmin=cmin,
                              cmax=cmax,
                              colorbar=dict(
                              #    title=cbar_label,
                                  x=cbar_loc,
                              )
                              )
        trace_z2 = go.Surface(z=list(z2_offset),
                              x=list(x_y),
                              y=list(y_y),
                              showscale=False,
                              surfacecolor=colorsurfz2,
                              colorscale=colorscale,
                              cmin=cmin,
                              cmax=cmax,
                              colorbar=dict(
                              #     title=cbar_label,
                                   x=cbar_loc,
                              )
                              )

        data = [trace_y1, trace_x1, trace_z1, trace_z2]

        layout = go.Layout(
                annotations=[
                    dict(text=r'$t='+str(round(time,isd))+unit+r'$',
                         x=textxy[0],
                         y=textxy[1],
                         showarrow=False
                        )],
                autosize=autosize,
                width=image_dim[0],
                height=image_dim[1],
                scene=dict(
                    camera=dict(
                        eye=dict(
                            x=viewpoint[0],
                            y=viewpoint[1],
                            z=viewpoint[2]
                               )),

                    xaxis=dict(
                        title=axestitle[0],
                        visible=visxyz[0]
                              ),
                    yaxis=dict(
                        title=axestitle[1],
                        visible=visxyz[1]
                              ),
                    zaxis=dict(
                        title=axestitle[2],
                        visible=visxyz[2]
                              )
                          ),
                margin=dict(
                            l=margin[0],
                            r=margin[1],
                            b=margin[2],
                            t=margin[3],
                           ),
                          )

        # plot the figure
        fig = go.Figure(data=data, layout=layout)
        # set filename of the figure
        filename = field+ '_{0:04d}'.format(itt) + "." + imageformat
        # display the figure html
        # po.plot(fig, image="png", auto_open=True,
        #          image_height=500, image_width=1000, filename=filename)
        #print the figure to file
        if not os.path.exists(figdir):
            os.mkdir(figdir)
        pio.write_image(fig, figdir + filename)


def plot_box(slice_obj,#slice_obj=pcn.read.slices()
             #acquire the slice objects
             fields=['uu1',], datadir='./data/', proc=-1, xyzplane=[],
             quiet=True, oldfile=False, 
             #select data to plot
             tstart=0., tend=1e38, islice=-1,
             #set image properties
             imageformat=".png", figdir='./images/',
             #set color parameters
             colorscale='Hot', norm='linear', 
             color_range=None, color_levels=None,
             #locate the box and axes
             viewpoint=(-1.35, -2.1, 0.5), offset=5., margin=(20,20,30,0),
             autosize=False, image_dim=(800,500),
             #handle axes properties
             visxyz=[True,True,True], axestitle=('x', 'y', 'z'), xyz=None,
             #color bar properties
             cbar_label = r'$u_x\,[{\rm km s}^{-1}]$', cbar_loc=1.,
             cbar_label_pos='top',
             #add text for time stamp                                     
             timestamp=False, textxy=(0,0), unit='', isd=2,  fontsize=25,
             ):

    #gd = pcn.read.grid(trim=True, quiet=True, datadir=datadir)
    ttmp=slice_obj.t[np.where(slice_obj.t<=tend)[0]]
    it=np.where(ttmp>=tstart)[0]
    if len(xyzplane)==0:
        for key in slice_obj.__dict__.keys():
            if key!='t':
                xyzplane.append(key)
    if len(xyzplane)<4:
        raise ValueError("xyzplane: rvid_box requires at least 4 surfaces.")
        
    if islice == -1:
        for itt in it:
            plot(
                 #field (or multiple todo) and 4 surfaces (extra todo) 
                 slice_obj, fields, xyzplane,
                 itt, it, quiet=quiet,  
                 #yz[i], xy[i], xz[i], xy2[i], 
                 #save output
                 figdir=figdir, imageformat=imageformat,
                 #set color parameters
                 norm=norm, colorscale=colorscale,
                 color_levels=color_levels,
                 color_range=color_range,
                 #locate the box and axes
                 viewpoint=viewpoint, offset=offset, margin=margin,
                 autosize=autosize, image_dim=image_dim,
                 #handle axes properties
                 visxyz=visxyz, axestitle=axestitle, xyz=xyz,
                 #color bar properties
                 cbar_label=cbar_label, cbar_loc=cbar_loc,
                 cbar_label_pos=cbar_label_pos,
                 #add text for time stamp                               
                 time=slice_obj.t[itt], textxy=textxy, unit=unit,
                 isd=isd, fontsize=25,  
                )
    else:
        plot(
             #field (or multiple todo) and 4 surfaces (extra todo) 
             slice_obj, fields, xyzplane,
             islice, [islice,], quiet=quiet,  
             #yz[i], xy[i], xz[i], xy2[i], 
             #save output
             figdir=figdir, imageformat=imageformat,
             #set color parameters
             norm=norm, colorscale=colorscale,
             color_levels=color_levels,
             color_range=color_range,
             #locate the box and axes
             viewpoint=viewpoint, offset=offset, margin=margin,
             autosize=autosize, image_dim=image_dim,
             #handle axes properties
             visxyz=visxyz, xyz=(gd.x, gd.y, gd.z), axestitle=axestitle,
             #color bar properties
             cbar_label=cbar_label, cbar_loc=cbar_loc,
             cbar_label_pos=cbar_label_pos,
             #add text for time stamp
             time=slice_obj.t[itt], textxy=textxy, unit=unit,
             isd=isd, fontsize=25,  
            )
