
def preparePlot(x_datasets, y_datasets,
                x_errors=False, y_errors=False,
                xmin=False, xmax=False, ymin=False, ymax=False,
                xlog=False, ylog=False,
                fig=False, ax=False,
                PROFILE=False):

  import numpy as np
  import matplotlib
  import matplotlib.pyplot as plt
  plt.ioff()
  from pencil import io
  from ...visu.internal import calc_lims
  from ...visu.internal import MinorSymLogLocator

  ## ESTIMATE LIMITS
  if np.any(x_errors) and np.shape(x_errors) == np.shape(x_datasets):
    xlim_min = calc_lims(val_min = np.min([np.min(array) for array in np.array(x_datasets)-np.array(x_errors)]))
    xlim_max = calc_lims(val_max = np.max([np.max(array) for array in np.array(x_datasets)+np.array(x_errors)]))
  else:
    xlim_min = calc_lims(val_min = np.min([np.min(array) for array in x_datasets]))
    xlim_max = calc_lims(val_max = np.max([np.max(array) for array in x_datasets]))

  if np.any(y_errors) and np.shape(y_errors) == np.shape(y_datasets):
    ylim_min = calc_lims(val_min = np.min([np.min(array) for array in np.array(y_datasets)-np.array(y_errors)]))
    ylim_max = calc_lims(val_max = np.max([np.max(array) for array in np.array(y_datasets)+np.array(y_errors)]))
  else:
    ylim_min = calc_lims(val_min = np.min([np.min(array) for array in y_datasets]))
    ylim_max = calc_lims(val_max = np.max([np.max(array) for array in y_datasets]))

  if ax:
    xlim_max = np.max([np.max(ax.get_xlim()), xlim_max]); xlim_min = np.min([np.min(ax.get_xlim()), xlim_min])
    ylim_max = np.max([np.max(ax.get_ylim()), ylim_max]); ylim_min = np.min([np.min(ax.get_ylim()), ylim_min])

  #pencil.io.debug_breakpoint()

  if not xmin: xmin = xlim_min
  if not xmax: xmax = xlim_max
  if not ymin: ymin = ylim_min
  if not ymax: ymax = ylim_max

  ## IF LOG: ENSHURE YOU ARE COVERING AT LEAST ONE ORDER OF MAGNITUDE
  #pencil.io.debug_breakpoint()
  while ylog and ymin > 0 and np.abs(np.int(np.floor(np.log10(ymin)))-np.int(np.floor(np.log10(ymax)))) < 1. and ymax > ymin and ymax*xmax != 0 :
    ymin = ymin - ymin * 0.01
    ymax = ymax + ymax * 0.01

  ## SWAP AXIS DIRECTION IF PURE NEGATIVE VALUE PLOT
  if xmin < 0 and xmax < 0 and xmin < xmax: tmp = xmax; xmax = xmin; xmin = tmp
  if ymin < 0 and ymax < 0 and ymin < ymax: tmp = ymax; ymax = ymin; ymin = tmp

  ## CREATE FIG AND AX IF NEEDED
  if (not fig) and (not ax):                 # if no fig/axis object is handed over, produce your own
    fig = plt.figure(figsize=PROFILE['fig_size'], facecolor=PROFILE['fig_background'])  # produce new figure if no axis is handed over
    ## SHOW PLOT?
    if matplotlib.get_backend().startswith('Tk'):
        print('## Showing plot..')
        plt.ion()
        plt.show()
    # Remove the plot frame lines. They are unnecessary chartjunk.
    ax = fig.add_subplot(111)

    ## prepare plot:
    ## remove axes lines if demanded
    ax.spines["top"].set_visible(PROFILE['ax_lines_top'])
    ax.spines["bottom"].set_visible(PROFILE['ax_lines_bottom'])
    ax.spines["right"].set_visible(PROFILE['ax_lines_right'])
    ax.spines["left"].set_visible(PROFILE['ax_lines_left'])

    if PROFILE['ax_ticks_top']: ax.get_xaxis().tick_top()
    if PROFILE['ax_ticks_bottom']: ax.get_xaxis().tick_bottom()
    if PROFILE['ax_ticks_left']: ax.get_yaxis().tick_left()
    if PROFILE['ax_ticks_right']: ax.get_yaxis().tick_right()

    ## set axes scales
    if xlog and ylog:                                   # Set x- and y-axis to log
      ax.set_xscale("log")
      ax.set_yscale("symlog", linthreshy=np.min([ymin,np.abs(ylim_min),np.abs(ylim_max),np.min(np.abs(y_datasets))])/10.)
    elif (not xlog) and ylog:
      ax.set_yscale("symlog", linthreshy=np.min([xmni,np.abs(ylim_min),np.abs(ylim_max),np.min(np.abs(y_datasets))])/10.)
    elif xlog and (not ylog):
      ax.set_xscale("log")
    else:
      print('!! ERROR: only log setting not defined so far! xlog='+str(xlog)+' and ylog='+str(ylog))  # nothing else implemented yet!

    ## set fontsize for ticks
    plt.setp(ax.get_xticklabels(), fontsize=PROFILE['ax_tick_fontsize'])
    plt.setp(ax.get_yticklabels(), fontsize=PROFILE['ax_tick_fontsize'])

    ## set ticks for symlog
    if ylog:
      yaxis = fig.gca().yaxis
      yaxis.set_minor_locator(MinorSymLogLocator(1e-1))

  ## ESTABLISH CALCULATED LIMITS FOR PLOT
  #pencil.io.debug_breakpoint()
  ax.set_xlim(left=xmin)
  ax.set_xlim(right=xmax)
  ax.set_ylim(bottom=ymin)
  ax.set_ylim(top=ymax)

  if ax == False:
      print("!! ERROR: Created axis object is False! Starting debug_breakpoint")
      io.debug_breakpoint()

  return (fig, ax)
