'''General purpose plotting for models.'''

import numpy as np
import matplotlib.pyplot as plt
import math


TIME_VAR = "time"


class DataPlotter(object):
  """
  Provides plots for data sources.
  The data source must provide the DataProvider interface.
    getVariable(variable_name) - returns np.array
      (Raises ValueError if the name is not found.)
  Many methods in DataPlotter take as argument list-of-names.
  """

  def __init__(self, data_provider):
    """
    :param Object data_provider: implements the DataProvider interface.
    """
    self._provider = data_provider

  def lines(self, names, yrange=None):
    """
    Does a stacked line plot with a common time axis.
    """
    num_plots = len(names)
    xvals = self._provider.getVariable(TIME_VAR)
    xmin = min(xvals)
    xmax = max(xvals)
    for idx in range(num_plots):
      plt.subplot(num_plots, 1, idx+1)
      frame = plt.gca()
      if idx < num_plots - 1:
        frame.axes.get_xaxis().set_visible(False)
      else:
        plt.xlabel(TIME_VAR)
      plt.ylabel(names[idx])
      yvals = self._provider.getVariable(names[idx])
      if yrange is None:
        ymin = min(yvals)*0.95
        ymax = max(yvals)*1.05
      else:
        ymin = yrange[0]
        ymax = yrange[1]
      locs,labels = yticks()
      yticks(locs, map(lambda x: "%.1f" % x, locs*1e9))
      plt.axis([xmin, xmax, ymin, ymax])
      plt.plot(xvals, yvals)
