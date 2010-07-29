# $Id: animate_interactive.py,v 1.5 2010-07-02 18:54:58 iomsn Exp $
#
# read 3D array and assemble an animation
#
# Author: Simon Candelaresi (iomsn@physto.se). 
# 
#

import numpy as np
import pylab as plt
import time
import os, sys # for making the movie
#from matplotlib.widgets import Slider, Button # for GUI
from matplotlib.backends.backend_tkagg import *
import mtTkinter as tk # for thread save GUI (from http://tkinter.unpythonic.net/wiki/mtTkinter)
import sys
import thread # for GUI

def animate_interactive(data, t = [], dimOrder = (0,1,2),
                        fps = 10.0, title = '', xlabel = '', ylabel = '',
                        cBar = 0, interpol = 'nearest', colorTable = 'hot',
                        aspectRatio = 'auto',
                        movieFile = '', keepImages = False):
    """
    read 3D array and assemble an animation.

    Options:

     data --- array with time dependent slices (3d)
     t --- time steps corresponding to the slices
     dimOrder --- ordering of the dimensions in the array (t,x,y)
     fps --- frames per second of the animation
     title --- title of the plot
     xlabel --- label of the x-axis
     ylabel --- label of the y-axis
     cBar --- determines how the colorbar should be change
              (0 - no cahnge; 1 - keep extreme values constant; 2 - change extreme values)
     interpol --- interpolation of the data points
     colorTable --- colortable for ht plot
     aspectRatio --- aspect ration of the plot
     movieFile --- the movie file where the animation should be saved to,
                   leave empty for no movie, requires 'mencoder' to be installed
     keepImages --- if True the images for the movie creation are not deleted
    """
    
    global tStep, sliderTime, pause
    
    
    # plot the current frame
    def plotFrame():
        global tStep, sliderTime
    
        ax.set_title('t={0}'.format(t[tStep]))
        
        if (cBar == 0):
            image.set_data(data[tStep,:,:])
        if (cBar == 1):
            image.set_data(data[tStep,:,:])
            colorbar.set_clim(vmin=data[tStep,:,:].min(), vmax=data[tStep,:,:].max())
        if (cBar == 2):
            image.set_data(data[tStep,:,:])
            colorbar.set_clim(vmin=data[tStep,:,:].min(), vmax=data[tStep,:,:].max())
            colorbar.update_bruteforce(data[tStep,:,:])

        #manager.canvas.draw()
        canvas.draw()
        #print 'draw', tStep
                
    
    # play the movie
    def play(threadName):               
        global tStep, sliderTime, pause, rootWindow
        
        pause = False
        while (tStep < nT) & (pause == False):        
            # write the image files for the movie
            if movieFile:
                plotFrame()
                frameName = movieFile + '%06d.png'%tStep
                fig.savefig(frameName)
                movieFiles.append(frameName)
            else:
                start = time.clock()
                # time slider
                sliderTime.set_val(t[tStep])
                # wait for the next frame (fps)
                while (time.clock() - start < 1.0/fps):
                    pass # do nothing                                        
            tStep += 1        
        tStep -= 1
            
    
    # call the play function as a seperate thread (for GUI)
    def play_thread():
        try:
            thread.start_new_thread(play, ("playThread", ))
        except:
            print "Error: unable to start play thread"


    def pausing():               
        global pause
        
        pause = True        


    def reverse():
        global tStep, sliderTime
        
        tStep -= 1
        if tStep < 0:
            tStep = 0
        #plotFrame()
        # plot the frame and update the time slider
        sliderTime.set_val(t[tStep])

        
    def forward():
        global tStep, sliderTime
        
        tStep += 1
        if tStep > len(t)-1:
            tStep = len(t)-1
        #plotFrame()
        # plot the frame and update the time slider
        sliderTime.set_val(t[tStep])

    
    pause = False
    
    # of no movie file is specified interacrive mode is switched on
    interactive = True
    if movieFile:
        interactive = False
        
    # check if the data has the right dimensions
    if (len(data.shape) != 3):
        print 'error: data dimensions are invalid: {0} instead of 3'.format(len(data.shape))
        return -1
    
    # transpose the data according to dimOrder
    unOrdered = data
    data = np.transpose(unOrdered, dimOrder)
    unOrdered = []
    
    # check if time array has the right length
    nT = len(t)
    if (nT != len(data[:,0,0])):
        print 'error: length of time array doesn\'t match length of data array'
        return -1
    
    # check if fps is positive
    if (fps < 0.0):
        print 'error: fps is not positive, fps = {0}'.format(fps)
        return -1

    # determine the size of the array
    nX = len(data[0,:,0])
    nY = len(data[0,0,:])
    
    # determine the minimum and maximum values of the data set
    dataMin = np.min(data)
    dataMax = np.max(data)
    
    # setup the plot and the frame for the GUI
    rootWindow = tk.Tk()
    fig = plt.Figure()
    ax = fig.add_subplot(111)
    
    # create the time slider
    if interactive == True:
        fig.subplots_adjust(bottom=0.5)
        #fig.axes.append(plt.axes([0.2, 0.1, 0.6, 0.03], axisbg='lightgoldenrodyellow'))
        #sliderTimeAxes = plt.axes([0.2, 0.1, 0.6, 0.03], axisbg='lightgoldenrodyellow')
        sliderTimeAxes = fig.add_axes([0.2, 0.1, 0.6, 0.03], axisbg='lightgoldenrodyellow')
        #sliderTimeAxes = fig.gca()
        #sliderTimeAxes = fig.axes[-1]
        sliderTime = plt.Slider(sliderTimeAxes, 'time', t[0], t[-1], valinit = 0.0)
        sliderTimeAxes.set_figure(fig)
        def update(val):
            global tStep
            # find the closest time step to the slider time value
            for i in range(len(t)):
                if t[i] < sliderTime.val:
                    tStep = i
            if (tStep != len(t)-1):
                if (t[tStep+1] - sliderTime.val) < (sliderTime.val - t[tStep]):
                    tStep += 1
            plotFrame()
        sliderTime.on_changed(update)
    
    canvas = FigureCanvasTkAgg(fig, master=rootWindow)
    tkwidget = canvas.get_tk_widget()
    tkwidget.pack(side=tk.TOP)
    
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    plane = np.zeros((nX,nY))

    # calibrate the displayed colors for the data range
    image = ax.imshow(plane, vmin=dataMin, vmax=dataMax, interpolation=interpol, cmap=colorTable, aspect=aspectRatio)
    colorbar = fig.colorbar(image)

    ## for real-time image display
    #manager = plt.get_current_fig_manager()
    #manager.show()
    canvas.show()


    tStep = 0
    if movieFile:
        movieFiles = []
        # start the animation
        play('noThread')

        # write the movie file        
        mencodeCommand = "mencoder 'mf://"+movieFile+"*.png' -mf type=png:fps="+np.str(fps)+" -ovc lavc -lavcopts vcodec=mpeg4:vhq:vbitrate=1800 -ffourcc MP4S -oac copy -o "+movieFile+".mpg"
        os.system(mencodeCommand)
        # clean up the image files
        if (keepImages == False):
            print 'cleaning up files'
            for fname in movieFiles:
                os.remove(fname)

    else:
        # set up the gui        
        plt.ion()
        # new GUI
        frame = tk.Frame(rootWindow)
        frame.pack(side=tk.BOTTOM)
        
        buttonPlay = tk.Button(frame, text="play", command=play_thread)
        buttonPlay.pack(side=tk.LEFT)
        buttonPause = tk.Button(frame, text="pause", command=pausing)
        buttonPause.pack(side=tk.LEFT)
        
        buttonReverse = tk.Button(frame, text="reverse", command=reverse)
        buttonReverse.pack(side=tk.LEFT)
        buttonForward = tk.Button(frame, text="forward", command=forward)
        buttonForward.pack(side=tk.LEFT)
                
        rootWindow.mainloop()
    print 'done'
