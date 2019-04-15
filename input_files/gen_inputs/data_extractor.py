from netCDF4 import Dataset
import numpy as np
from gl_dataset import *
from pyproj import Proj
import matplotlib.pyplot as plt
from scipy import interpolate
from pyproj import Proj
from scipy.spatial import cKDTree
import fenics as fc

font = {'size'   : 17}
plt.rc('font', **font)

class DataExtractor():

    def __init__(self):

        ### Load the precip. and temp. data
        ########################################################################

        # Precip. data
        precip_data = Dataset('/home/jake/flow_model/paleo_data/modern/Box_precip.nc')
        # Temp. data
        temp_data = Dataset('/home/jake/flow_model/paleo_data/modern/Box_temp.nc')
        # Lat, lon
        self.lats = precip_data.variables['lat'][:]
        self.lons = precip_data.variables['lon'][:]
        # Reference elevation
        self.dem = precip_data.variables['dem'][:]
        # Temperature
        self.temp = temp_data.variables['Temperature'][:]
        # Precipitation
        self.precip = precip_data.variables['SMB'][:]

        
        ### Compute 30 year average temperature across each month as modern
        ########################################################################
        self.monthly_temp_modern = []
        self.monthly_precip_modern = []

        for i in range(12):
            self.monthly_temp_modern.append(np.average(self.temp[140:170,i], axis = 0))
            self.monthly_precip_modern.append(np.average(self.precip[140:170,i], axis = 0))

        # KD Tree for finding lat / lon points in data set near flowline points
        self.tree = cKDTree(np.c_[self.lons.ravel(),  self.lats.ravel()])
        # Standard map coord. projection
        self.projection = Proj(init='epsg:3413')
        # Bed data
        self.bed_data = GLDataset('bed')
        # Surface data
        self.surface_data = GLDataset('surface')


    """
    Saves the modern day temp and precip along the flowline to files.
    """
    def write_mesh(self, file_name, flow_lats, flow_lons, res = 1000.0):

        xs, ys = self.projection(self.lons.flatten(), self.lats.flatten())
        xs_flow, ys_flow = self.projection(flow_lons, flow_lats)

        bed_xs = self.bed_data.bed_xarray
        bed_ys = self.bed_data.bed_yarray


        ### Compute evenly spaced points along flowline at given resolution
        ########################################################################

        # Convert lat, lon coords to map coords
        xs, ys = self.projection(flow_lons, flow_lats)
        # Create a function that maps from distance to x, y position
        segment_lens = np.sqrt((xs[1:] - xs[:-1])**2 + (ys[1:] - ys[:-1])**2)
        Ls = np.insert(segment_lens.cumsum(), 0, 0.)
        x_interp = interpolate.interp1d(Ls, xs)
        y_interp = interpolate.interp1d(Ls, ys)
        # Get x, y coordinates every res meters
        Ls_res = np.arange(Ls[0], Ls[-1], res)
        xs_res = x_interp(Ls_res)
        ys_res = y_interp(Ls_res)
        # Recompute lat, lon coords at the given resolution
        flow_lons, flow_lats = self.projection(xs_res, ys_res, inverse=True)


        ### Create an hdf5 file and mesh
        ########################################################################

        mesh = fc.IntervalMesh(len(Ls_res)-1, 0., 1.)
        V = fc.FunctionSpace(mesh, "CG", 1)
        V_r = fc.FunctionSpace(mesh, "R", 0)
        V_dg = fc.FunctionSpace(mesh, 'DG', 0)
        # Bed fenics function
        bed_func = fc.Function(V)
        # Surface fenics function
        surface_func = fc.Function(V)
        # Precip fenics function
        precip_func = fc.Function(V)
        # Temp fenics function
        temp_func = fc.Function(V)
        # Domain length fenics function
        domain_len_func = fc.Function(V_r)
        # Create HDF5 output file
        out_file = fc.HDF5File(mesh.mpi_comm(), file_name, "w")
        # Write mesh to output file
        out_file.write(mesh, "mesh")


        ### Get bed data along flowline
        ########################################################################

        bed_res = 3500.
        bed_Ls = np.linspace(Ls[0], Ls[-1], Ls[-1] / bed_res, endpoint = True)
        xs_bed = x_interp(bed_Ls)
        ys_bed = y_interp(bed_Ls)
        bed_interp = interpolate.UnivariateSpline(bed_Ls, self.bed_data.getInterpolatedValue(xs_bed, ys_bed))
        flow_bed = bed_interp(Ls_res)
        # Write bed data
        bed_func.vector()[:] = np.ascontiguousarray(flow_bed[::-1])
        out_file.write(bed_func, "B")
        # Write domain length
        domain_len_func.vector()[:] = Ls_res[-1]
        out_file.write(domain_len_func, "domain_len")


        ### Get the surface data along the flowline
        ########################################################################

        surface_interp = interpolate.UnivariateSpline(bed_Ls, self.surface_data.getInterpolatedValue(xs_bed, ys_bed))
        flow_surface = surface_interp(Ls_res)
        # Write surface data
        surface_func.vector()[:] = np.ascontiguousarray(flow_surface[::-1])
        out_file.write(surface_func, "S")
        # Write domain length
        fc.plot(surface_func)
        fc.plot(bed_func)
        plt.show()
        

        ### Calculate temp and precip along flowline
        ########################################################################

        # Get indexes of closest points to each flowline point
        dd, ii = self.tree.query(np.c_[flow_lons,  flow_lats], k = 50)
        # Unique indexes
        ii = np.unique(ii)
        # Convert indexes to row, col format
        row_indexes = ii / self.lons.shape[1]
        col_indexes = np.mod(ii, self.lons.shape[1])
        # Get lat, lon coords and indexes
        lats = self.lats[row_indexes, col_indexes]
        lons = self.lons[row_indexes, col_indexes]

        # Interpolate DEM over region around flowline
        dem_interp = interpolate.LinearNDInterpolator(np.c_[lons, lats], self.dem[row_indexes, col_indexes])
        flow_dem = dem_interp(np.c_[flow_lons, flow_lats])

        S_ref = fc.Function(V)
        S_ref.vector()[:] = np.ascontiguousarray(flow_dem[::-1])
        out_file.write(S_ref, "S_ref")

        for i in range(12):
            # Get modern temp and precip for given month
            temp = self.monthly_temp_modern[i]
            precip = self.monthly_precip_modern[i]

            # Interpolate temp over region around flowline
            temp_interp = interpolate.LinearNDInterpolator(np.c_[lons, lats], temp[row_indexes, col_indexes])
            # Interpolate precip over region around flowline
            precip_interp = interpolate.LinearNDInterpolator(np.c_[lons, lats], precip[row_indexes, col_indexes])
            # Get interpolated temp and precip
            flow_temp = temp_interp(np.c_[flow_lons, flow_lats])
            flow_precip = precip_interp(np.c_[flow_lons, flow_lats])
            # Convert precip to meters per year
            flow_precip = (flow_precip / 1000.) * 12.

            precip_func.vector()[:] = np.ascontiguousarray(flow_precip)
            temp_func.vector()[:] = np.ascontiguousarray(flow_temp)

            out_file.write(temp_func, "T" + str(i))
            out_file.write(precip_func, "P" + str(i))
            fc.plot(precip_func)

        plt.show()


        ### Write initial ice sheet configuration data
        ########################################################################
        out_file.close()

