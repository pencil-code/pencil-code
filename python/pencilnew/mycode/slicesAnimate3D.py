#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 26 15:04:37 2018

@author: root
"""
"""
animate the contour of a field in three surfaces of the box.
input:
    *field*:    the field you want to animate
    *proc*:     the processor number that slices data from
    *filename*: the name of the generated video
output:
    a video: which shows the animation of slices, 
             saved in the current directory
             
# =============================================================================
# call signature: slicesAnimate3D(field='bb1', proc=0, filename='bb1_xy3D.mp4')                                                   
# =============================================================================
    
# =============================================================================
# attention: close the animation figure to end the function
# =============================================================================
"""
def slicesAnimate3D(*args, **kwargs):
    animate_temp = __SlicesAnimate3D__(*args, **kwargs)
    animate_temp.__read__()
    animate_temp.__animate__()
    return animate_temp

class __SlicesAnimate3D__(object):
    """
    self arguments:
    *t*:          time series
    *slicesXY*: the 2D array series of xy surface
    *slicesXZ*: the 2D array series of xz surface
    *slicesYZ*: the 2D array series of yz surface
    *slicesXY2*: the 2D array series of xy2 surface
    *filename*:   used to store the output animate video
    """
    def __init__(self, field='bb1', proc=-1, filename=None):
        
        self.field = field
        self.proc = proc
        if(filename == None):
            self.filename = field+'3D.mp4'
    
    def __read__(self):
        from ..read import slices, grid
        self.slices = slices(field=self.field,extension=['xy','xz','yz','xy2'], proc=self.proc)
        self.grid = grid(proc=self.proc,trim=True)
        
        
    def __animate__(self):
        import numpy as np
        import matplotlib.pyplot as plt
        import matplotlib.animation as animation
        
        
    
        fig = plt.figure()
        ax = fig.gca(projection='3d')
        #title = plt.title()
#        X,Y = np.meshgrid(self.grid.x, self.grid.y)
#        Z = getattr(getattr(self.slices, 'xy'), self.field)
#        ax.contourf(X, Y, Z[0,:,:], zdir='z', offset = max(self.grid.z), cmap='plasma')
#        
#        X,Z = np.meshgrid(self.grid.x, self.grid.z)
#        Y = getattr(getattr(self.slices, 'xz'), self.field)
#        ax.contourf(X, Y[0,:,:], Z, zdir='y', offset = min(self.grid.y), cmap='plasma')
#        
#        Y,Z = np.meshgrid(self.grid.y, self.grid.z)
#        X = getattr(getattr(self.slices, 'yz'), self.field)
#        ax.contourf(X[0,:,:], Y, Z, zdir='x', offset = max(self.grid.x), cmap='plasma')
#        
#        X,Y = np.meshgrid(self.grid.x, self.grid.y)
#        Z = getattr(getattr(self.slices, 'xy2'), self.field)
#        ax.contourf(X, Y, Z[0,:,:], zdir='z', offset = min(self.grid.z)-0.7*self.grid.Lz,
#                    cmap='plasma')
#        
#        # Create cubic bounding box to simulate equal aspect ratio
#        max_range = np.array([self.grid.Lx, self.grid.Ly, 1.7*self.grid.Lz]).max()
#        Xb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][0].flatten() + 0.5*(self.grid.x.max()+self.grid.x.min())
#        Yb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][1].flatten() + 0.5*(self.grid.y.max()+self.grid.y.min())
#        Zb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][2].flatten() + 0.5*(self.grid.x.max()+self.grid.x.min()-0.7*self.grid.Lz)
#        # Comment or uncomment following both lines to test the fake bounding box:
#        for xb, yb, zb in zip(Xb, Yb, Zb):
#           ax.plot([xb], [yb], [zb], 'w')
#        
#        def init():
#            frame.set_data([self.slice_data[0,:,:]])
#            #title
#            return frame,
#        
        def updatefig(i):
            ax.clear()
            X,Y = np.meshgrid(self.grid.x, self.grid.y)
            Z = getattr(getattr(self.slices, 'xy'), self.field)
            ax.contourf(X, Y, Z[i,:,:], zdir='z', offset = max(self.grid.z), cmap='plasma')
            
            X,Z = np.meshgrid(self.grid.x, self.grid.z)
            Y = getattr(getattr(self.slices, 'xz'), self.field)
            ax.contourf(X, Y[i,:,:], Z, zdir='y', offset = min(self.grid.y), cmap='plasma')
            
            Y,Z = np.meshgrid(self.grid.y, self.grid.z)
            X = getattr(getattr(self.slices, 'yz'), self.field)
            ax.contourf(X[i,:,:], Y, Z, zdir='x', offset = max(self.grid.x), cmap='plasma')
            
            X,Y = np.meshgrid(self.grid.x, self.grid.y)
            Z = getattr(getattr(self.slices, 'xy2'), self.field)
            ax.contourf(X, Y, Z[i,:,:], zdir='z', offset = min(self.grid.z)-0.7*self.grid.Lz,
                        cmap='plasma')
            
            # Create cubic bounding box to simulate equal aspect ratio
            max_range = np.array([self.grid.Lx, self.grid.Ly, 1.7*self.grid.Lz]).max()
            Xb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][0].flatten() + 0.5*(self.grid.x.max()+self.grid.x.min())
            Yb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][1].flatten() + 0.5*(self.grid.y.max()+self.grid.y.min())
            Zb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][2].flatten() + 0.5*(self.grid.x.max()+self.grid.x.min()-0.7*self.grid.Lz)
            for xb, yb, zb in zip(Xb, Yb, Zb):
                ax.plot([xb], [yb], [zb], 'w')
        
            ax.set_title('t = %f' % self.slices.t[i] )
            return ax
        
        ani = animation.FuncAnimation(fig, updatefig,len(self.slices.t), interval=50, blit=False)

        ani.save(self.filename,fps=5,dpi=400, writer='ffmpeg')
#   
        #ax.set_zlim(min(self.grid.z)-3*self.grid.Lz, max(self.grid.z))
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        
        plt.show()
