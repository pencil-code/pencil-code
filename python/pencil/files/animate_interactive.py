# $Id: animate_interactive.py iomsn Exp $
#
# Assemble a 2D animation from a 3D array.
#
# Author: Simon Candelaresi (iomsn@physto.se, iomsn1@googlemail.com).
# 
#

import numpy as np
import pylab as plt
import time
import os # for making the movie
from matplotlib.backends.backend_tkagg import *
import thread # for GUI

def animate_interactive(data, t = [], dimOrder = (0,1,2),
                        fps = 10.0, title = '', xlabel = 'x', ylabel = 'y',
                        fontsize = 24,
                        cBar = 0,
                        #interpol = 'nearest', colorTable = 'hot', aspectRatio = 'auto',
                        arrowsX = np.array(0), arrowsY = np.array(0), arrowsRes = 10,
                        arrowsPivot = 'mid', arrowsWidth = 0.002, arrowsScale = 5,
                        arrowsColor = 'black', plotArrowsGrid = False,
                        movieFile = '', bitrate = 1800, keepImages = False,
                        **kwimshow):
    """
    Assemble a 2D animation from a 3D array.

    call signature::
    
      animate_interactive(data, t = [], dimOrder = (0,1,2),
                        fps = 10.0, title = '', xlabel = 'x', ylabel = 'y',
                        fontsize = 24,
                        cBar = 0,
                        #interpol = 'nearest', colorTable = 'hot', aspectRatio = 'auto',
                        arrowsX = np.array(0), arrowsY = np.array(0), arrowsRes = 10,
                        arrowsPivot = 'mid', arrowsWidth = 0.002, arrowsScale = 5,
                        arrowsColor = 'black', plotArrowsGrid = False,
                        movieFile = '', keepImages = False,
                        **kwimshow)
    
    Assemble a 2D animation from a 3D array. *data* has to be a 3D array who's
    time index has the same dimension as *t*. The time index of *data* as well
    as its x and y indices can be changed via *dimOrder*.
    
    Keyword arguments:
    
      *dimOrder*: [ (i,j,k) ]
        Ordering of the dimensions in the data array (t,x,y).
        
     *fps*:
       Frames per second of the animation.
       
     *title*:
       Title of the plot.
       
     *xlabel*:
       Label of the x-axis.
       
     *ylabel*:
       Label of the y-axis.
       
     *fontsize*:
       Font size of the title, x and y label.
       
     *cBar*: [ 0 | 1 | 2 ]
       Determines how the colorbar changes:
       (0 - no cahnge; 1 - keep extreme values constant; 2 - change extreme values).
              
     *arrowsX*:
       Data containing the x-component of the arrows.
       
     *arrowsy*:
       Data containing the y-component of the arrows.
       
     *arrowsRes*:
       Plot every arrowRes arrow.
       
     *arrowsPivot*: [ 'tail' | 'middle' | 'tip' ]
       The part of the arrow that is at the grid point; the arrow rotates
       about this point.
       
     *arrowsWidth*:
       Width of the arrows.
       
     *arrowsScale*:
       Scaling of the arrows.
       
     *arrowsColor*:
       Color of the arrows.
       
     *plotArrowsGrid*: [ False | True ]
       If 'True' the grid where the arrows are aligned to is shown.
     
     *movieFile*: [ None | string ]
       The movie file where the animation should be saved to.
       If 'None' no movie file is written. Requires 'mencoder' to be installed.
     
     *bitrate*:
       Bitrate of the movie file. Set to higher value for higher quality.
       
     *keepImages*: [ False | True ]
       If 'True' the images for the movie creation are not deleted.
     
     **kwimshow:
       Remaining arguments are identical to those of pylab.imshow. Refer to that help.
    """


    global tStep, sliderTime, pause
    
    
    # plot the current frame
    def plotFrame():
        global tStep, sliderTime
        
        if movieFile:
            ax.set_title(title+r'$\quad$'+r'$t={0}$'.format(t[tStep]), fontsize = fontsize)
        
        if (cBar == 0):
            image.set_data(data[tStep,:,:])
        if (cBar == 1):
            image.set_data(data[tStep,:,:])
            colorbar.set_clim(vmin=data[tStep,:,:].min(), vmax=data[tStep,:,:].max())
        if (cBar == 2):
            image.set_data(data[tStep,:,:])
            colorbar.set_clim(vmin=data[tStep,:,:].min(), vmax=data[tStep,:,:].max())
            colorbar.update_bruteforce(data[tStep,:,:])
            
        if plotArrows:
            arrows.set_UVC(U = arrowsX[tStep,::arrowsRes,::arrowsRes], V = arrowsY[tStep,::arrowsRes,::arrowsRes])

        #manager.canvas.draw()
                
    
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
    def play_thread(event):
        try:
            thread.start_new_thread(play, ("playThread", ))
        except:
            print "Error: unable to start play thread"


    def pausing(event):               
        global pause
        
        pause = True        


    def reverse(event):
        global tStep, sliderTime
        
        tStep -= 1
        if tStep < 0:
            tStep = 0
        # plot the frame and update the time slider
        sliderTime.set_val(t[tStep])

        
    def forward(event):
        global tStep, sliderTime
        
        tStep += 1
        if tStep > len(t)-1:
            tStep = len(t)-1
        # plot the frame and update the time slider
        sliderTime.set_val(t[tStep])

    
    pause = False
    plotArrows = False
    
    # check if the data has the right dimensions
    if (len(data.shape) != 3):
        print 'error: data dimensions are invalid: {0} instead of 3'.format(len(data.shape))
        return -1
        
    # transpose the data according to dimOrder
    unOrdered = data
    data = np.transpose(unOrdered, dimOrder)
    unOrdered = []
    
    # check if arrows should be plotted
    if len(arrowsX.shape) == 3:
        # transpose the data according to dimOrder
        unOrdered = arrowsX
        arrowsX = np.transpose(unOrdered, dimOrder)
        unOrdered = []
        if len(arrowsY.shape) == 3:
            # transpose the data according to dimOrder
            unOrdered = arrowsY
            arrowsY = np.transpose(unOrdered, dimOrder)
            unOrdered = []
            
            # check if the dimensions of the arrow arrays match each other
            if ((len(arrowsX[:,0,0]) != len(arrowsY[:,0,0])) or (len(arrowsX[0,:,0]) != len(arrowsY[0,:,0])) or (len(arrowsX[0,0,:]) != len(arrowsY[0,0,:]))):
                print 'error: dimensions of arrowX do not match with dimensions of arrowY'
                return -1
            else:
                plotArrows = True
    
    # check if time array has the right length
    nT = len(t)
    if (nT != len(data[:,0,0])):
        print 'error: length of time array doesn\'t match length of data array'
        return -1
        if plotArrows:
            if (nT != len(arrowX[:,0,0]) or nT != len(arrowX[:,0,0])):
                print 'error: length of time array doesn\'t match length of arrows array'
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
    
    # setup the plot
    if movieFile:
        width = 8
        height = 7
        plt.rc("figure.subplot", left=(90/72.27)/width)
        plt.rc("figure.subplot", right=(width-20/72.27)/width)
        plt.rc("figure.subplot", bottom=(60/72.27)/height)
        plt.rc("figure.subplot", top=(height-20/72.27)/height)
        fig = plt.figure(figsize=(width, height))        
        ax = plt.axes([0.11, 0.1, .90, .85])
    else:
        width = 8
        height = 7
        plt.rc("figure.subplot", left=(90/72.27)/width)
        plt.rc("figure.subplot", right=(width-20/72.27)/width)
        plt.rc("figure.subplot", bottom=(60/72.27)/height)
        plt.rc("figure.subplot", top=(height-20/72.27)/height)
        fig = plt.figure(figsize=(width, height))
        ax = plt.axes([0.1, 0.25, .85, .70])
        #ax = fig.add_subplot(111)
    
    ax.set_title(title, fontsize = fontsize)
    ax.set_xlabel(xlabel, fontsize = fontsize)
    ax.set_ylabel(ylabel, fontsize = fontsize)
    plane = np.zeros((nX,nY))

    # calibrate the displayed colors for the data range
    image = ax.imshow(plane, vmin=dataMin, vmax=dataMax, **kwimshow)
    colorbar = fig.colorbar(image)
    
    # plot the arrows
    # TODO: add some more options
    if plotArrows:
        # prepare the mash grid where the arrows will be drawn
        arrowGridX, arrowGridY = np.meshgrid(np.arange(0, len(data[0,:,0]), arrowsRes), np.arange(0, len(data[0,0,:]), arrowsRes))
        arrows = ax.quiver(arrowGridX, arrowGridY, arrowsX[0,::arrowsRes,::arrowsRes], arrowsY[0,::arrowsRes,::arrowsRes], units = 'width', pivot = arrowsPivot, width = arrowsWidth, scale = arrowsScale, color = arrowsColor)
        # plot the grid for the arrows
        if plotArrowsGrid == True:
            ax.plot(arrowGridX, arrowGridY, 'k.')
        

    # for real-time image display
    #manager = plt.get_current_fig_manager()
    #manager.show()


    tStep = 0
    if movieFile:
        movieFiles = []
        # start the animation
        play('noThread')

        # write the movie file        
        mencodeCommand = "mencoder 'mf://"+movieFile+"*.png' -mf type=png:fps="+np.str(fps)+" -ovc lavc -lavcopts vcodec=mpeg4:vhq:vbitrate="+np.str(bitrate)+" -ffourcc MP4S -oac copy -o "+movieFile+".mpg"
        os.system(mencodeCommand)
        # clean up the image files
        if (keepImages == False):
            print 'cleaning up files'
            for fname in movieFiles:
                os.remove(fname)

    else:
        # set up the gui        
        plt.ion()

        axPlay = plt.axes([0.1, 0.05, 0.15, 0.05], axisbg='lightgoldenrodyellow')
        buttonPlay = plt.Button(axPlay, 'play', color='lightgoldenrodyellow', hovercolor='0.975')
        buttonPlay.on_clicked(play_thread)
        axPause = plt.axes([0.3, 0.05, 0.15, 0.05], axisbg='lightgoldenrodyellow')
        buttonPause = plt.Button(axPause, 'pause', color='lightgoldenrodyellow', hovercolor='0.975')
        buttonPause.on_clicked(pausing)
    
        axReverse = plt.axes([0.5, 0.05, 0.15, 0.05], axisbg='lightgoldenrodyellow')
        buttonReverse = plt.Button(axReverse, 'reverse', color='lightgoldenrodyellow', hovercolor='0.975')
        buttonReverse.on_clicked(reverse)
        axForward = plt.axes([0.7, 0.05, 0.15, 0.05], axisbg='lightgoldenrodyellow')
        buttonForward = plt.Button(axForward, 'forward', color='lightgoldenrodyellow', hovercolor='0.975')
        buttonForward.on_clicked(forward)
        
        # create the time slider
        fig.subplots_adjust(bottom=0.2)
        sliderTimeAxes = plt.axes([0.2, 0.12, 0.6, 0.03], axisbg='lightgoldenrodyellow')
        sliderTime = plt.Slider(sliderTimeAxes, 'time', t[0], t[-1], valinit = 0.0)
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
    
        plt.show()
    print 'done'
