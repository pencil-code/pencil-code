#!/usr/bin/env python

#
# $Id: multi_slices.py,v 1.1 2009-08-27 09:57:23 iomsn Exp $
# Show specified slices over time interactively.
#

#
# important information for developers:
# This program uses two different indexations: one runs over all buttons,
# the other one runs over the subplots (axes). For index conversion there are
# the two button2axes and axes2button. In a following version i will try to
# get rid of this.
# The following are twi lists eplaining which lists and functions are using
# button or axes indexing.
#
# lists:
# gui.buttonFlag    button
# gui.buttonObj     button
# plot.graphOptions button
# plot.axes         axes
# plot.images       axes
# 
# functions:
# button2axes       button
# axes2button       axes
# newPlot           axes
# replot            axes
# buttonPressed     button

import numpy
import pylab as P
import pencil_old as pc
import sys
import time
import Tkinter as Tk # for GUI
import Pmw # for GUI

# some constants
BTNS_X = 6; BTNS_Y = 4  # number of displayed buttons in x- and y-direction

# class containing global variables for the GUI
class guiClass:
    def __init__(self):
        root = []   # root window with the console
        
        # 0 = not pressed; 1 = pressed; 2 = currently selected
        self.buttonFlag = list(0 for i in range(0,BTNS_X*BTNS_Y,1)) # idx: button
        self.buttonPre = -1  # previously pressed button
        self.buttonSelected = -1  # currently selected button
        self.buttonFrame = []    # frame for the button widget
        self.buttonObj = []  # the actual button object, idx: button
        
        # sub class which contains the option menues
        class optionMenuClass:
            def __init__(self):
                self.quantity = '<none>'
                self.component = '<none>'
                self.slices = '<none>'
        self.optionMenu = optionMenuClass()

        # widgets inside the gui window
        self.widgetOptions = []
gui = guiClass()

# class containing global variables for the plots
class plotClass:
    def __init__(self):
        # sub class which contains the graph options
        class graphOptionsClass:
            def __init__(self):
                self.quantity = '<none>'
                self.component = '<none>'
                self.slices = '<none>'
        self.graphOptions = []  # idx: button
        for each in range(BTNS_X*BTNS_Y): self.graphOptions.append(graphOptionsClass())
        
        self.figure = 0     # main figure containing the subplots (axes)
        self.axes = []      # axes with the subplots, idx: axes
        # contains the image displayed in the axes
        self.images = []    # idx: axes
        self.axesX = 0      # number of axes in x direction
        self.axesY = 0      # number of axes in y direction
        self.axWidth = 0    # width of the axes
        self.axHeight = 0   # height of the axes
        
        self.sliderAxes = []    # axes for the time slider
        self.slider = []        # the time slider object
plot = plotClass()

# class containing global variables for the data
class dataClass:
    def __init__(self):
        self.timeStep = 0    # the current time step
        self.slices = {}     # contains the data from the slices
        self.t = []          # all time steps
        self.dim = {}        # system dimansion
        self.param = {}      # system parameters
        self.loaded = set([])    # set containing informaion about which data has already been loaded
        self.possibleData = []   # list containing all possible field+extension combinations
        self.null = [[0,1],[2,3]]   # contains filling data
data = dataClass()


################################################################################
# data stuff
################################################################################

# load the data from the slice files
def loadData(field, extension):
    global data
    
    # check if data has already been loaded
    if (field+extension in data.loaded) == False:
        data.slices[field+extension], data.t = pc.read_slices(datadir = 'data', field = field, extension = extension, proc = -1)
        data.dim[field+extension] = pc.read_dim() # system dimensions
        data.param[field+extension] = pc.read_param(quiet=True) # system parameters
        data.loaded.add(field+extension)
    
# converts 'quantity' and 'component' into 'field'
def convertField(quantity, component):
    field = quantity    
    if (component == 'x'):
        field += '1'
    elif (component == 'y'):
        field += '2'
    elif (component == 'z'):
        field += '3'
    return field

# check if field data exists
def dataExists(field, extension):
    global data
    
    if field+extension in data.possibleData:
        return True
    else: return False


################################################################################
# plot stuff (pylab)
################################################################################

# convert button index into axes index
def button2axes(button):
    global gui, plot

    # find out the x- and y-coordinates of the buttons
    # NB: here i use integer division
    buttonY = int(button/BTNS_X)
    buttonX = button - buttonY*BTNS_X
    axes = buttonX + buttonY*plot.axesX
    return axes
    
# convert axes index into buton index
def axes2button(axes):
    global gui, plot

    # find out the x- and y-coordinates of the axes
    # NB: here i use integer division
    axesY = int(axes/plot.axesX)
    axesX = axes - axesY*plot.axesX
    button = axesX + axesY*BTNS_X
    return button
    
# determine the geometry of the subplots (axes)
def geometry():
    global gui, plot
    
    plot.axesX = 0
    plot.axesY = 0
    for btnX in range(BTNS_X):
        for btnY in range(BTNS_Y):
            if (gui.buttonFlag[btnX+btnY*BTNS_X] != 0):
                if (btnX > plot.axesX): plot.axesX = btnX
                if (btnY > plot.axesY): plot.axesY = btnY
    plot.axesX += 1
    plot.axesY += 1
    plot.axWidth = 0.8/plot.axesX
    plot.axHeigth = 0.8/plot.axesY


 # make time slider if necessary
def makeSlider():
    global data, plot
    
    if plot.sliderAxes == []:
        plot.sliderAxes = P.axes([0.2, 0.05, 0.6, 0.03], axisbg='lightgoldenrodyellow') # make place for a slider
    if plot.slider == []:
        plot.slider = P.Slider(plot.sliderAxes, 'time', data.t[0], data.t[-1], valinit=data.t[0]) # make a slider for the time
        plot.slider.on_changed(updatePlots)   # define the function to be executed when the slider value is changed
        P.draw()


# restructures the axes in the main figure
# only needed when axes are added or removed
def restructure():
    global gui, plot    
    
    # determine the geometry of the subplots (axes)
    geometry()

    # create the main figure
    if (plot.figure == 0):
        plot.figure = P.figure()
    
    # delete first all previous axes to avoid overdrawing
    if (plot.axes != []):
        for i in range(len(plot.axes)): plot.figure.delaxes(plot.axes[i])            
    plot.axes = []
    plot.images = []
    
    # create sub plots (axes)
    for ay in range(plot.axesY):
        for ax in range(plot.axesX):
            pos = [0.1+ax*plot.axWidth*1.1, 1.0-(ay+1)*plot.axHeigth*1.1, plot.axWidth, plot.axHeigth]
            plot.axes.append(plot.figure.add_axes(pos))
            axes = ax+ay*plot.axesX
                        
            # put an initial image in the axes
            button = axes2button(axes)
            field = convertField(plot.graphOptions[button].quantity, plot.graphOptions[button].component)
            # append entry into images list
            plot.images.append(P.imshow(data.null))
            if dataExists(field, plot.graphOptions[button].slices):
                newPlot(axes)
                makeSlider()
    P.draw()
        

# creates a new plot in the axes 'axes'
# NB: be aware of the different indexations (button <-> axes)
def newPlot(axes):
    global data, plot
    
    button = axes2button(axes)
    field = convertField(plot.graphOptions[button].quantity, plot.graphOptions[button].component)
    plot.images[axes] = plot.axes[axes].imshow(data.slices[field+plot.graphOptions[button].slices][data.timeStep,...], cmap=P.cm.hot, aspect='auto', origin='lower') # plot the slice
    P.draw()


# replot one specific axes
def replot(axes):
    global data, plot
    
    # since axes and images use a different indexing to plot.graphOptions
    # a conversion is needed
    #ti = time.time()
    button = axes2button(axes)
    #print "replot axes2button(axes): dt = ", time.time()-ti
    
    #ti = time.time()
    field = convertField(plot.graphOptions[button].quantity, plot.graphOptions[button].component)
    #print "replot convertField(): dt = ", time.time()-ti
    
    # change the image of this axes
    #ti = time.time()
    plot.images[axes].set_data(data.slices[field+plot.graphOptions[button].slices][data.timeStep])
    #print "replot plot.images[axes].set_data(): dt = ", time.time()-ti


# update the plots when the slider is used
def updatePlots(val):
    global data, plot
    
    print("updatePlot")
    #ti = time.time()
    data.timeStep = round((plot.slider.val-data.t[0])/(data.t[-1]-data.t[0])*(len(data.t)-1))
    #print "updatePlots round(): dt = ", time.time()-ti
    # replot all axes
    for axes in range(plot.axesX*plot.axesY):
        button = axes2button(axes)
        replot(axes)
    ti = time.time()
    P.draw()
    #print "updatePlots P.draw: dt = ", time.time()-ti
    
################################################################################
# console stuff (Tkinter)
################################################################################

def buttonPressed(btn):
    global gui    

    restruct = 0
    if gui.buttonFlag[btn] == 0:
        gui.buttonObj[btn].config(relief=Tk.SUNKEN, background='red', activebackground='red')
        gui.buttonFlag[btn] = 2
        restruct = 1
    elif gui.buttonFlag[btn] == 1:
        gui.buttonObj[btn].config(background='red', activebackground='red')
        gui.buttonFlag[btn] = 2
    elif gui.buttonFlag[btn] == 2:
        gui.buttonObj[btn].config(relief=Tk.RAISED, background='#e0dfde', activebackground='#e0dfde')
        gui.buttonFlag[btn] = 0
        restruct = 1
    if (gui.buttonPre > -1) and (gui.buttonPre != btn):
        gui.buttonObj[gui.buttonPre].config(background='#e0dfde', activebackground='#e0dfde')
        if gui.buttonFlag[gui.buttonPre] == 2:
            gui.buttonFlag[gui.buttonPre] = 1

    gui.buttonSelected = btn
    gui.buttonPre = btn

    # display options for selected graph in the option menu
    gui.optionMenu.quantity.setvalue(plot.graphOptions[btn].quantity)
    gui.optionMenu.component.setvalue(plot.graphOptions[btn].component)
    gui.optionMenu.slices.setvalue(plot.graphOptions[btn].slices)
    
    if (restruct == 1):
        restructure()   # restructure the axes


def changeQuantity(void):
    global gui, plot

    if (gui.buttonSelected != -1):
        axes = button2axes(gui.buttonSelected)
        plot.graphOptions[gui.buttonSelected].quantity = gui.optionMenu.quantity.getcurselection()
        field = convertField(gui.optionMenu.quantity.getcurselection(), gui.optionMenu.component.getcurselection())
        # check if field data exists
        if dataExists(field, gui.optionMenu.slices.getcurselection()):
            # load the data
            loadData(field, gui.optionMenu.slices.getcurselection())
            # mke time slider if necessary
            makeSlider()
            # assign the data to the axes
            newPlot(axes)
    

def changeComponent(void):
    global gui, plot

    if (gui.buttonSelected != -1):
        axes = button2axes(gui.buttonSelected)
        plot.graphOptions[gui.buttonSelected].component = gui.optionMenu.component.getcurselection()
        field = convertField(gui.optionMenu.quantity.getcurselection(), gui.optionMenu.component.getcurselection())
        # check if field data exists
        if dataExists(field, gui.optionMenu.slices.getcurselection()):
            # load the data
            loadData(field, gui.optionMenu.slices.getcurselection())
            # mke time slider if necessary
            makeSlider()
            # assign the data to the axes
            newPlot(axes)


def changeSlices(void):
    global gui, plot

    if (gui.buttonSelected != -1):
        axes = button2axes(gui.buttonSelected)
        plot.graphOptions[gui.buttonSelected].slices = gui.optionMenu.slices.getcurselection()
        field = convertField(gui.optionMenu.quantity.getcurselection(), gui.optionMenu.component.getcurselection())
        # check if field data exists
        if dataExists(field, gui.optionMenu.slices.getcurselection()):
            # load the data
            loadData(field, gui.optionMenu.slices.getcurselection())
            # mke time slider if necessary
            makeSlider()
            # assign the data to the axes
            newPlot(axes)


################################################################################
# main program
################################################################################

# initialize window with the console
gui.root = Tk.Tk()
gui.root.title('pvid_slices.py')


# initialize set which contains all possible field+extension combinations
# if the pencil code gets addition possible quantities for the slices
# this set must be changed accordingly
possibleQuantities = ['bb1', 'bb2', 'bb3', 'aa1', 'aa2', 'aa3', 'uu1', 'uu2', 'uu3', 'oo1', 'oo2', 'oo3', 'b2', 'u2', 'o2', 'divu', 'lnrho', 'ss', 'shock', 'lnTT', 'ionization', 'Qrad', 'Isurf', 'lncc', 'ecr']
for i in range(len(possibleQuantities)):
    data.possibleData.append(possibleQuantities[i]+'xy')
    data.possibleData.append(possibleQuantities[i]+'yz')
    data.possibleData.append(possibleQuantities[i]+'xz')
    data.possibleData.append(possibleQuantities[i]+'xy2')
del(possibleQuantities)


# buttons
for btnY in range(0,BTNS_Y,1):
    gui.buttonFrame.append(Tk.Frame(gui.root))
    gui.buttonFrame[btnY].pack(side=Tk.TOP)
    for btnX in range(0,BTNS_X,1):
        #q1 = lambda: buttonPressed(btnX, btnY)
        gui.buttonObj.append(Tk.Button(gui.buttonFrame[btnY],command=lambda btn=(btnX+btnY*BTNS_X): buttonPressed(btn), relief=Tk.RAISED))
        gui.buttonObj[btnX+btnY*BTNS_X].pack(side=Tk.LEFT)


# drop down menues widget
gui.widgetOptions = Pmw.Group(tag_text='physical quantities and slices')
gui.widgetOptions.pack(side = 'bottom', padx = 5, pady = 3)

gui.optionMenu.quantity = Pmw.OptionMenu(gui.widgetOptions.interior(),
    labelpos = 'w',
    label_text = 'quantity:',
    items = ['bb', 'b2', 'aa', 'uu', 'u2', 'oo', 'o2', 'divu', 'lnrho', 'ss', 'shock', 'lnTT', 'ionization', 'Qrad', 'Isurf', 'lncc', 'ecr', '<none>'],
    command = changeQuantity,
    menubutton_width = 8,
)
gui.optionMenu.quantity.pack(side = 'left', padx = 5, pady = 3)
gui.optionMenu.quantity.invoke('<none>')

gui.optionMenu.component = Pmw.OptionMenu(gui.widgetOptions.interior(),
    labelpos = 'w',
    label_text = 'component:',
    items = ['x', 'y', 'z', '<none>'],
    command = changeComponent,
    menubutton_width = 8,
)
gui.optionMenu.component.pack(side = 'left', padx = 5, pady = 3)
gui.optionMenu.component.invoke('<none>')

gui.optionMenu.slices = Pmw.OptionMenu(gui.widgetOptions.interior(),
    labelpos = 'w',
    label_text = 'slice:',
    items = ['xy', 'xz', 'yz', 'xy2', '<none>'],
    command = changeSlices,
    menubutton_width = 8,
)
gui.optionMenu.slices.pack(side = 'left', padx = 5, pady = 3)
gui.optionMenu.slices.invoke('<none>')
