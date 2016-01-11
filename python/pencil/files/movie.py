# $Id: movie.py 11852 2009-11-09 07:39:23Z rplasson $
#
# write movie from data
#
# Author: R. Plasson (rplasson@nordita.org). 
# 
#

def make_movie(data,tmin=0,tmax=0,dt=1,output="/tmp/out.avi",dim=-1,min=-1,max=1,stride=1):
    """
    read 2D slices from data and assemble an animation in a movie.

    Quickly written from the example at http://matplotlib.sourceforge.net/faq/howto_faq.html
 
    Options:

     data -- function, data(t) returns an array to plot
     tmin -- starting time
     tmax -- end time
     dt -- step time
     output -- output filename
     min -- minimum data value
     max -- maximum data value
     def --  stride for 3D plot
    """
    import os
    import pylab as P
    import tempfile as T
    import numpy as N
    from mpl_toolkits.mplot3d import Axes3D

    
    files=""
    flist=[]
    fig = P.figure()
    grid=data(0).shape
    if dim==-1:  # autodetect dimension
        dim=len(grid)
    if dim==3:   # prepare mesh in 3D case
        x,y=N.meshgrid(N.arange(grid[0]),N.arange(grid[1]))
        ax= Axes3D(fig)
    else:
        ax = fig.add_subplot(111)
        
    for t in range(tmin,tmax,dt):
        ax.cla()
        if dim==1:
            ax.plot(data(t),)
            P.ylim([min,max])
        elif dim==3:
            ax.plot_surface(x,y,data(t),rstride=stride, cstride=stride, cmap=P.cm.jet)
            ax.set_zlim3d([min,max])
        else:
            ax.imshow(data(t),vmin=min,vmax=max)
        tmp=T.NamedTemporaryFile(suffix=".png")
        fname=tmp.name
        flist+=[tmp]
        #print 'Saving frame', t # Python 2
        print('Saving frame '+t)
        fig.savefig(fname)
        files+=fname+","
            
    #print 'Making movie animation.mpg - this make take a while' # Python 2
    print('Making movie animation.mpg - this make take a while')
    os.system("mencoder mf:/"+files[:-1]+" -mf type=png:fps=24 -ovc lavc -lavcopts vcodec=wmv2 -oac copy -o "+output)
    for i in flist:
        del i            # Purge temp files
