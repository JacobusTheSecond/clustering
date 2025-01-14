from abc import ABC, abstractmethod
import numpy as np
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
import csv
import os
import math

import klcluster as kl

class DriftersSolver(ABC):
    def __init__(self, drifterFiles):
        self.clustercurves = None
        self.curves_per_cluster = []
        self.datacurves = []
        for filepath in drifterFiles:
            with open(filepath, "r") as f:
                curveData = np.array(list(csv.reader(f, delimiter=" "))).astype(float)
                if len(curveData) > 1:
                    self.datacurves.append(curveData)
        self.rossbyRadii = {}
        with open(os.path.dirname(os.path.realpath(__file__))+'/../data_drifters/rossrad.dat', 'r') as file:
            for line in file:
                line = line.strip().split()
                self.rossbyRadii[(float(line[0]), float(line[1]))] = float(line[3])
                #print(line)
    @abstractmethod
    def solve():
        pass
    
    def plotInput(self):
        self.__plotCurves([self.datacurves], ["blue"], [0.4])

    def plotResults(self):
        if (not self.clustercurves):
            raise Exception("No result found, consider calling solve() before")
        self.__plotCurves([self.clustercurves], ["red"], [0.7])

    def plotInputAndResult(self):
        self.__plotCurves([self.datacurves, self.clustercurves], ["blue", "red"], [0.4, 0.8])

    def plotCurve(self):
        #self.curves_per_cluster.sort(reverse=True, key=lambda x: len(x[0]))
        for i in range(2):
            lat_min, lon_min, lat_max, lon_max = 180, 180, -180, -180
            for curve in self.curves_per_cluster[i][0]:
                for c in curve:
                
                    lat, lon = self.worldToLL(c[0], c[1], c[2])
                    lat_min = min(lat, lat_min)
                    lon_min = min(lon, lon_min)
                    lat_max = max(lat, lat_max)
                    lon_max = max(lon, lon_max)
            self.__plotCurves([self.datacurves, self.curves_per_cluster[i][0], [self.curves_per_cluster[i][1]]], ["gray", "blue", "red"], [0.2, 0.4, 0.9], [lon_min-2, lon_max+2, lat_min-2, lat_max+2])

    def __plotCurves(self, curve_list, color_list=["blue"], size_list = [0.5], custom_window=None):
        ax = plt.axes(projection=ccrs.PlateCarree())
        
        if custom_window is None:
            arrow_scale = 1
            ax.set_global()
        else:
            arrow_scale = (custom_window[1]-custom_window[0])/180
            ax.set_extent(custom_window)
        ax.coastlines()
        #ax.stock_img()
        plt.tight_layout()
        ax.get_figure().canvas.manager.set_window_title("Drifters")

        def add_arrow(line, size=0.00001, color=None, arrow_scale=arrow_scale):
            """
            add an arrow to a line.

            line:       Line2D object
            position:   x-position of the arrow. If None, mean of xdata is taken
            direction:  'left' or 'right'
            size:       size of the arrow in fontsize points
            color:      if None, line color is taken.
            """
            if color is None:
                color = line.get_color()

            xdata = line.get_xdata()
            ydata = line.get_ydata()

            if len(xdata) < 2:
                return

            index = len(xdata) - 1
            x = xdata[index]
            y = ydata[index]
            dx = (xdata[index]-xdata[index-1])*0.01
            dy = (ydata[index]-ydata[index-1])*0.01
            norm = (dx**2 + dy**2)**0.5
            dx_norm = (dx / norm) * 0.1 
            dy_norm = (dy / norm) * 0.1

            arrow_scale *= line.get_linewidth()
            # Add the arrow to the line
            plt.arrow(x, y, dx_norm, dy_norm, facecolor=color, edgecolor=color,
                    head_width=arrow_scale, head_length=arrow_scale, zorder=5 if color=="red" else 1)
            #plt.arrow(x, y, dx, dy, facecolor=color, edgecolor=color,
            #head_width=1, head_length=1, transform=ccrs.Geodetic())

            # line.axes.annotate('',
            #     [x + dx, y + dy],
            #     [x, y],
            #     arrowprops=dict(width=0, headwidth=3, headlength=3, color=color),
            #     annotation_clip=True
            # )
        for curves, color, width in zip(curve_list, color_list, size_list):
            
            for i, curve in enumerate(curves):
                if curve.size < 3:
                    continue
                lat, lon = self.worldToLL(curve[:,0], curve[:,1], curve[:,2])

                line = plt.plot(lon, lat,
                        color=color, linewidth=width,
                        transform=ccrs.Geodetic())[0]
                
                add_arrow(line)

                if i % 100 == 0:
                    print(f"{i}/{len(curves)}")
                i += 1
                # if i > 600:
                #     break


        plt.show()

    def llToWorld(self, lat,lon):
        lat, lon = np.deg2rad(lat), np.deg2rad(lon)
        R = 6371000 # radius of the earth
        x = R * np.cos(lat) * np.cos(lon)
        y = R * np.cos(lat) * np.sin(lon)
        z = R *np.sin(lat)
        return np.array([x,y,z])

    def worldToLL(self, x, y, z):
        R = 6371000 # radius of the earth
        lat = np.arcsin(z / R) / np.pi * 180
        lon = np.arctan2(y, x) / np.pi * 180
        return np.array([lat, lon])
    def getRossbyRadius(self, lat, lon):
        lon = (lon + 360) % 360
        lat = math.floor(lat) + 0.5
        lon = math.floor(lon) + 0.5
        return self.rossbyRadii.get((lat, lon), 1)
    
    def drawInputField(self):
        self.__drawVectorField(self.datacurves, "blue")

    def drawResultField(self):
        #interpolated_coords = [self.clustercurves[0]]  # Start with the first coordinate
        threshold = 100000
        new_curves = []
        for clustercurve in self.clustercurves:
            interpolated_curve = []
            for i in range(1, len(clustercurve)):
                x1, y1, z1 = clustercurve[i - 1]
                x2, y2, z2 = clustercurve[i]
                
                # Compute the Euclidean distance
                distance = np.sqrt((x2 - x1)**2 + (y2 - y1)**2 + (z2 - z1)**2)
                
                if distance > threshold:
                    # Number of points to interpolate
                    num_points = int(np.ceil(distance / threshold))
                    for t in np.linspace(0, 1, num_points + 1)[1:]:  # Avoid repeating the start point
                        interpolated_curve.append((
                            x1 + t * (x2 - x1),
                            y1 + t * (y2 - y1),
                            z1 + t * (z2 - z1)
                        ))
                else:
                    interpolated_curve.append((x2, y2, z2))
            new_curves.append(interpolated_curve)

        self.__drawVectorField(new_curves, "red")

    def __drawVectorField(self, curves, color):
        # # Sample trajectory data with only lat/lon points
        # trajectory_data = [
        #     [(-74, 40), (-30, 50), (2, 48)],   # NYC to Paris
        #     [(-0.1, 51.5), (60, 55), (139, 35)],  # London to Tokyo
        #     [(144.9, -37.8), (120, -10), (80, 0)],  # Melbourne to Southeast Asia
        # ]

        # Initialize lists to hold points and vectors

        lon, lat, u, v = [], [], [], []

        # Calculate vectors based on differences between consecutive points
        for path in curves:
            for i in range(len(path) - 1):
                # Get start and end points for vector calculation
                x, y, z = path[i]
                lat1, lon1 = self.worldToLL(x, y, z)
                x, y, z = path[i+1]

                lat2, lon2 = self.worldToLL(x, y, z)
            
                # Calculate the vector components as differences in lon/lat
                delta_lon = lon2 - lon1
                delta_lat = lat2 - lat1
                #rescale = math.sqrt(delta_lat**2+delta_lon**2)
                #delta_lat /= rescale
                #delta_lon /= rescale
                #print(delta_lat, delta_lon)
            
                # Store the midpoint for plotting the vector
                mid_lon = (lon1 + lon2) / 2
                mid_lat = (lat1 + lat2) / 2
            
                # Append midpoint and vector to lists
                lon.append(mid_lon)
                lat.append(mid_lat)
                u.append(delta_lon)
                v.append(delta_lat)

        # Convert lists to arrays
        lon = np.array(lon)
        lat = np.array(lat)
        u = np.array(u)
        v = np.array(v)

        # Define the grid parameters
        grid_size = 5  # Grid cell size in degrees
        lon_bins = np.arange(-180, 180, grid_size)
        lat_bins = np.arange(-90, 90, grid_size)
        lon_centers = lon_bins[:-1] + grid_size / 2
        lat_centers = lat_bins[:-1] + grid_size / 2

        # Initialize arrays to hold the averaged vector components
        u_avg = np.zeros((len(lat_centers), len(lon_centers)))
        v_avg = np.zeros((len(lat_centers), len(lon_centers)))
        counts = np.zeros((len(lat_centers), len(lon_centers)))  # Count of vectors in each cell

        # Aggregate vectors in each cell
        for i in range(len(lon)):
            lon_idx = np.digitize(lon[i], lon_bins) - 1
            lat_idx = np.digitize(lat[i], lat_bins) - 1
        
            # Only aggregate if indices are within the grid
            if 0 <= lon_idx < len(lon_centers) and 0 <= lat_idx < len(lat_centers):
                u_avg[lat_idx, lon_idx] += u[i]
                v_avg[lat_idx, lon_idx] += v[i]
                counts[lat_idx, lon_idx] += math.sqrt(u[i]**2+v[i]**2)

        # Average the vectors in each cell
        u_avg = np.divide(u_avg, counts, where=counts != 0)
        v_avg = np.divide(v_avg, counts, where=counts != 0)
        print(u_avg)

        # Mask zero vectors by setting them to NaN (so they won't be plotted)
        u_avg[counts == 0] = np.nan
        v_avg[counts == 0] = np.nan

        # Create a map with Cartopy
        fig, ax = plt.subplots(figsize=(12, 6), subplot_kw={'projection': ccrs.PlateCarree()})
        ax.set_global()
        ax.coastlines()
        ax.add_feature(cfeature.BORDERS, linestyle=':')

        # Add gridlines
        gl = ax.gridlines(draw_labels=True, linewidth=0.5, color='gray', alpha=0.5)
        gl.top_labels = gl.right_labels = False

        for lon_line in lon_bins:
            ax.plot([lon_line, lon_line], [lat_bins[0], lat_bins[-1]], transform=ccrs.PlateCarree(), color='gray', linewidth=0.5, linestyle='--')
        for lat_line in lat_bins:
            ax.plot([lon_bins[0], lon_bins[-1]], [lat_line, lat_line], transform=ccrs.PlateCarree(), color='gray', linewidth=0.5, linestyle='--')


        # Plot the averaged vector field, omitting zero vectors
        lon_grid, lat_grid = np.meshgrid(lon_centers, lat_centers)
        quiver = ax.quiver(
            lon_grid, lat_grid, u_avg, v_avg,
            transform=ccrs.PlateCarree(), color=color, headaxislength=3, headlength=3, scale=0.2, scale_units='xy', width=0.0012#, scale=20, scale_units='xy', width=0.005
        )
        #for path in self.datacurves:
        #    path_lon = [point[0] for point in path]
        #    path_lat = [point[1] for point in path]
        #    ax.plot(path_lon, path_lat, transform=ccrs.PlateCarree(), color='red', linestyle='--', marker='o', markersize=5)

        # Optionally, add a reference arrow for scale
        ax.quiverkey(quiver, X=0.9, Y=-0.1, U=1, label="1 unit vector", labelpos='E')

        # Set the title
        plt.title("Averaged Global Vector Field in Grid Cells (Zero Vectors Omitted)")

        # Display the plot
        plt.show()

    


class KlClusterDriftersSolver(DriftersSolver):
    def __init__(self, drifterFiles, rossbyRadius=False):
        super().__init__(drifterFiles)

        self.curves = kl.Curves()
        
        dc = []
        for curvedata in self.datacurves:
            lat, lon = self.worldToLL(curvedata[0][0], curvedata[0][1], curvedata[0][2])
            #if lat > 0 and lon > 0 and lon < 100:
            dc.append(curvedata)
        self.datacurves = dc
        for i, curvedata in enumerate(self.datacurves):
            #print(curvedata)
            curve = kl.Curve(curvedata, f"Curve_{i}")
            weights = curve.getWeights()
            if rossbyRadius:
                for j, p in enumerate(curve):
                    lat, lon = self.worldToLL(p[0], p[1], p[2])
                    weights[j] = (self.getRossbyRadius(lat, lon)/3)+10
                    #if weights[j] != 500:
                    #print(weights[j])
            curve.setWeights(weights)
            self.curves.add(curve)

        self.DELTA = 50000
        # self.freeDELTA = 300000
        # self.simpDELTA = 5000
        self.freeDELTA = 400000
        self.simpDELTA = 20000
        self.COMPLEXITY = 20
        self.ROUNDS = 1

        print(f"Inititalizing {len(self.curves)} curves")
        self.cc = kl.CurveClusterer()
        if rossbyRadius:
            self.cc.initCurvesDiffDelta(self.curves, self.simpDELTA, self.freeDELTA)
        else:
            #self.cc.initCurves(self.curves, self.DELTA)
            self.cc.initCurvesDiffDelta(self.curves, self.simpDELTA, self.freeDELTA)

        
        self.simplifiedCurves = self.cc.getSimplifications()

    def lowerbound(self,withShow=False):
        num = self.cc.greedyIndependent(self.COMPLEXITY,withShow)
        print(f"Size of Greedy Independent: {num}")
        return num

    def solutionSize(self,withShow = False):
        self.clusters = self.cc.greedyCover(self.COMPLEXITY, self.ROUNDS, withShow)

        filterCount = 0
        return len(self.clusters)

    def solve(self, onlyRelevantClusters = False, withShow = False):
        self.clusters = self.cc.greedyCover(self.COMPLEXITY, self.ROUNDS, withShow)

        filterCount = 0
        self.curves = self.cc.getCurves()
        # for i in range(len(self.curves)):
        #     print(self.curves[i][0][0], self.simplifiedCurves[i][0][0], "????")

        # convert cluster centers to np array
        self.clustercurves = []

        for i in range(len(self.curves)):
            print(self.curves[i][0][0], self.simplifiedCurves[i][0][0])



        self.curves_per_cluster = []

        for cluster in self.clusters:
            
            curve = []


            center = cluster.center()
            start = int(center.start.value)
            end = int(center.end.value)
            curve = self.simplifiedCurves[center.curve]

            if onlyRelevantClusters:
                if end-start <= 1 or len(cluster.values()) <= 1:
                    filterCount += 1
                    continue
            curvedata = []
            for i in range(start, end+1):
                if (i < 0 or i >= len(curve)):
                    raise Exception("Index out of bounds")
                # print(curve[i])                
                curvedata.append(curve[i].values)
            all_curves = []
            for c in cluster:
                start = int(c.start.value)
                end = int(c.end.value)
                curve = self.curves[c.curve]
                curvedata_ = []
                #print(curve)
                # for c1 in self.curves:
                #     print(len(c1))
                start = int(self.cc.mapToBase(c.curve, c.start).value)
                end = int(self.cc.mapToBase(c.curve, c.end).value)
                # print(start, end, len(curve), c.curve)
                # print(len([r.values for r in curve]))
                for i in range(start, end+1):
                    if (i < 0 or i >= len(curve)):
                        raise Exception("Index out of bounds")
                    # print(curve[i])                
                    curvedata_.append(curve[i].values)
                #all_curves.append(np.array(self.simplifiedCurves[c.curve]))
                all_curves.append(np.array(curvedata_))#self.simplifiedCurves[c.curve]))
            self.clustercurves.append(np.array(curvedata))

            self.curves_per_cluster.append((all_curves, np.array(curvedata)))

        self.curves = self.cc.getCurves()
        # for i in range(len(self.curves)):
        #     print(self.curves[i][0][0], self.simplifiedCurves[i][0][0], "????")
        # print(len(self.curves), len(self.simplifiedCurves))
        # print(len(self.curves[0]), len(self.simplifiedCurves[0]))
        if onlyRelevantClusters:
            print(f"Filtered out {filterCount} center curves")
