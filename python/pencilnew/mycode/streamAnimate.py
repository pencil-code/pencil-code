#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 27 10:09:25 2018

@author: root
"""
"""
animate the contour of a field in a surface.
input:
    *vectorfield*:     the field you want to animate which listed in video.in 
                 and must be a vector
    *extension*: the slice surface you want to animate
    *proc*:      the processor number that slices data from
    *filename*:  the name of the generated video
output:
    a video: which shows the animation of slices, 
             locates in the current directory
             
# =============================================================================
# call signature: streamAnimate(vectorfield='bb', extension='xy', proc=0, 
                                                filename='bb_xy_stream.mp4')
# =============================================================================
             
# =============================================================================
# attention: close the animation figure to get out the function
# =============================================================================
"""
def streamAnimate(*args, **kwargs):
    animate_temp = __StreamAnimate__(*args, **kwargs)
    animate_temp.__read__()
    animate_temp.__animate__()
    return animate_temp

class __StreamAnimate__(object):
    
    def __init__(self, vectorfield='bb', extension='xy', proc=-1, filename=None):
        
        self.vectorfield = vectorfield
        self.extension = extension
        self.proc = proc
        if(filename == None):
            self.filename = 'stream_'+vectorfield+'_'+extension+'.mp4'
        
    
    def __read__(self):
        from ..read import slices, grid
        
        self.grid = grid(proc=self.proc,trim=True)
        if(self.extension == 'xy'):
            self.slices = slices(field=[self.vectorfield+'1',self.vectorfield+'2'], 
                                 extension='xy',proc=self.proc)
            self.__U__ = getattr(self.slices.xy, self.vectorfield+'1')
            self.__V__ = getattr(self.slices.xy, self.vectorfield+'2')
            self.axislabel = ['X','Y']
            self.__X__ = self.grid.x
            self.__Y__ = self.grid.y
            
        elif(self.extension == 'xz'):
            self.slices = slices(field=[self.vectorfield+'1',self.vectorfield+'3'],
                                 extension='xz',proc=self.proc)
            self.__U__ = getattr(self.slices.xz, self.vectorfield+'1')
            self.__V__ = getattr(self.slices.xz, self.vectorfield+'3')
            self.axislabel = ['X','Z']
            self.__X__ = self.grid.x
            self.__Y__ = self.grid.z
            
        elif(self.extension == 'yz'):
            self.slices = slices(field=[self.vectorfield+'2',self.vectorfield+'3'],
                                 extension='yz',proc=self.proc)
            self.__U__ = getattr(self.slices.yz, self.vectorfield+'2')
            self.__V__ = getattr(self.slices.yz, self.vectorfield+'3')
            self.axislabel = ['Y','Z']
            self.__X__ = self.grid.y
            self.__Y__ = self.grid.z
            
        elif(self.extension == 'xy2'):
            self.slices = slices(field=[self.vectorfield+'1',self.vectorfield+'2'],
                                 extension='xy2',proc=self.proc)
            self.__U__ = getattr(self.slices.xy, self.vectorfield+'1')
            self.__V__ = getattr(self.slices.xy, self.vectorfield+'2')
            self.axislabel = ['X','Y']
            self.__X__ = self.grid.x
            self.__Y__ = self.grid.y
            
        else:
            print('wrong extension, please check the input argument.')
            return 0
        
    def __animate__(self):
        import numpy as np
        import matplotlib.pyplot as plt
        import matplotlib.animation as animation
        
        fig = plt.figure()
        ax = fig.gca()
        ax.set_aspect(1)
        strength = np.sqrt(self.__U__*self.__U__+self.__V__*self.__V__)
        lw = 5*strength/strength.max()
        X,Y = np.meshgrid(self.__X__,self.__Y__)
#        def init():
#            frame.set_data([self.slice_data[0,:,:]])
#            #title
#            return frame,
#        
        def updatefig(i):
            ax.clear()
            ax.streamplot(X,Y,self.__U__[i,:,:],self.__V__[i,:,:],linewidth=lw[i,:,:])
            ax.set_xlim(X.min(),X.max())
            ax.set_ylim(Y.min(),Y.max())
            #frame.set_data(self.slice_data[i,:,:])
            ax.set_title('t = %f' % self.slices.t[i] )
            return ax
        
        ani = animation.FuncAnimation(fig, updatefig,len(self.slices.t), interval=50, blit=False)
        #ani = animation.FuncAnimation(fig, updatefig,30, interval=50, blit=False)


        ani.save(self.filename,fps=5,dpi=60, writer='ffmpeg')

        plt.show()
