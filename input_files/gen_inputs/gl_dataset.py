from pylab import sqrt, linspace
from scipy.interpolate import RectBivariateSpline
import numpy as np
import h5py

'''
Class: Dataset
Argument list: name of dataset, pen type(used for plotting)
Purpose: This is the class of datasets. This will store velocity, smb, etc. This takes the Velocity in X and Y direction
and makes one dataset of just Velocity. This velocity dataset ONLY stores the magnitude but not direction.

Dependencies: pylabs sqrt and linspace, RectBivariateSplint, numpy
Creator: James Stauder
Date created:2/23/18
Last edited: 3/2/18
'''

class GLDataset:
    def __init__(self, name):
        self.name = name

        proj_x0 = -637925
        proj_x1 = 864625
        proj_y0 = -657675
        proj_y1 = -3349425
        x1 = 1670
        y1 = 2991
        bed_xarray = linspace(proj_x0, proj_x1, x1, endpoint=True)
        bed_yarray = linspace(proj_y1, proj_y0, y1, endpoint=True)

        if self.name == 'velocity':
            self.data, self.vx, self.vy = self.setData(name)

            self.vxInterp = RectBivariateSpline(bed_xarray, bed_yarray, np.flipud(self.vx).transpose())
            self.vyInterp = RectBivariateSpline(bed_xarray, bed_yarray, np.flipud(self.vy).transpose())

        else:
            self.data = self.setData(name)

        self.bed_xarray = bed_xarray
        self.bed_yarray = bed_yarray

        self.interp = RectBivariateSpline(bed_xarray, bed_yarray, np.flipud(self.data).transpose())


    def setData(self, name):
        dataFile = h5py.File('/media/drive1/refactoring/data/GreenlandInBedCoord.h5', 'r')
        if name == 'velocity':
            vx = dataFile['VX'][:]
            vy = dataFile['VY'][:]
            data = sqrt(vx ** 2 + vy ** 2)
            dataFile.close()
            return data, vx, vy
        else:
            data = dataFile[name][:]
            dataFile.close()
            return data

    def getInterpolatedValue(self, xPosition, yPosition):
        return self.interp(xPosition, yPosition, grid = False)
