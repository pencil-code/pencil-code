# -*- coding: utf-8 -*-
"""
Created on Wed Apr 25 16:53:41 2018

@author: root
"""
"""
animate the contour of a field in a surface.
input:
    *field*:     the field you want to animate
    *extension*: the slice surface you want to animate
    *proc*:      the processor number that slices data from
    *filename*:  the name of the generated video
output:
    a video: which shows the animation of slices, 
             locates in the current directory

# =============================================================================
# call signature: slicesAnimate(field='bb1', extension='xy', proc=0, 
                                                        filename='bb1_xy.mp4')
# =============================================================================
             
# =============================================================================
# attention: close the animation figure to get out the function
# =============================================================================
"""
def slicesAnimate(*args, **kwargs):
    animate_temp = __SlicesAnimate__(*args, **kwargs)
    animate_temp.__read__()
    animate_temp.__animate__()
    return animate_temp

class __SlicesAnimate__(object):
    
    def __init__(self, field='bb1', extension='xy', proc=-1, filename=None):
        
        self.field = field
        self.extension = extension
        self.proc = proc
        if(filename == None):
            self.filename = field+'_'+extension+'.mp4'
    
    def __read__(self):
        from ..read import slices, grid
        self.slices = slices(field=self.field,extension=self.extension, proc=self.proc)
        self.grid = grid(proc=self.proc,trim=True)    # 坐标的数值分割
        
        
    def __animate__(self):
        import numpy as np
        import matplotlib.pyplot as plt
        import matplotlib.animation as animation
    
        fig = plt.figure()
        ax = fig.gca()
        #title = plt.title()
        #ax.contourf(self.slice_data[0,:,:],cmap='plasma')
        X,Y = np.meshgrid(self.grid.x, self.grid.y)
        Z = getattr(getattr(self.slices, self.extension), self.field)
#        def init():
#            frame.set_data([self.slice_data[0,:,:]])
#            #title
#            return frame,
#        
        def updatefig(i):
            ax.clear()
            ax.contourf(X, Y, Z[i,:,:],cmap='plasma')
            ax.set_aspect(1)
            #frame.set_data(self.slice_data[i,:,:])
            ax.set_title('t = %f' % self.slices.t[i] )
            return ax
        
        ani = animation.FuncAnimation(fig, updatefig,len(self.slices.t), interval=50, blit=False)

        ani.save(self.filename,fps=5,dpi=80, writer='ffmpeg')

        plt.show()
