from abc import ABC, abstractmethod
import numpy as np
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
import csv
import os
from frechet import euclidean, FastDiscreteFrechetMatrix
import math
import similaritymeasures
# from shapely import frechet_distance
# from frechetdist import frdist
import timeit
import ssl
ssl._create_default_https_context = ssl._create_unverified_context
import klcluster as kl

class DriftersSolver(ABC):
    def __init__(self, drifterFiles, points=10000):
        points = int(points)
        self.points = points
        self.clustercurves = None
        self.relevantClustercurves = None
        self.curves_per_cluster = []
        self.relevantCurves_per_cluster = []
        self.datacurves = []
        self.iterations = 1
        self.execution_time = 0
        self.init_time = 0
        
        for filepath in drifterFiles:
            with open(filepath, "r") as f:
                curveData = np.array(list(csv.reader(f, delimiter=" "))).astype(float)
                if len(curveData) > 1:
                    print(points)
                    if points > len(curveData):
                        self.datacurves.append(curveData)
                        #print(len(curveData), curveData)
                    else:
                        self.datacurves.append(curveData[0:len(curveData)-points])
                        break
                    points -= len(curveData)

                    #print(points)
        # self.rossbyRadii = {}
        # with open(os.path.dirname(os.path.realpath(__file__))+'/../data_drifters/rossrad.dat', 'r') as file:
        #     for line in file:
        #         line = line.strip().split()
        #         self.rossbyRadii[(float(line[0]), float(line[1]))] = float(line[3])
        #         #print(line)
        #self.frechet = FastDiscreteFrechetMatrix(euclidean)
    @abstractmethod
    def solve():
        pass
    
    def plotInput(self, filename="Input"):
        self.__plotCurves([self.datacurves], ["blue"], [0.4])

    def plotResults(self, relevant=False,num=-1, filename="Results"):
        if (not self.clustercurves):
            raise Exception("No result found, consider calling solve() before")
        filename += "Relevant" if relevant else "All"
        self.__plotCurves([self.relevantClustercurves if relevant else (self.clustercurves if num == -1 else self.clustercurves[0:num])], ["red"], [0.7], filename=filename)

    def plotInputAndResult(self, relevant=False, filename="InputAndResults"):
        filename += "Relevant" if relevant else "All"
        self.__plotCurves([self.datacurves, self.relevantClustercurves if relevant else self.clustercurves], ["blue", "red"], [0.2, 0.5])

    def plotCurve(self, relevant=False, filename="SingleCurve"):
        curves_per_cluster = self.relevantCurves_per_cluster if relevant else self.curves_per_cluster
        for i in range(3):
            lat_min, lon_min, lat_max, lon_max = 180, 180, -180, -180
            for curve in curves_per_cluster[i][0]:
                for c in curve:
                
                    lat, lon = self.worldToLL(c[0], c[1], c[2])
                    lat_min = min(lat, lat_min)
                    lon_min = min(lon, lon_min)
                    lat_max = max(lat, lat_max)
                    lon_max = max(lon, lon_max)
            self.__plotCurves([self.datacurves, curves_per_cluster[i][0], [curves_per_cluster[i][1]]], ["gray", "blue", "red"], [0.2, 0.4, 0.9], [lon_min-2, lon_max+2, lat_min-2, lat_max+2], filename+str(i+1))

    def __plotCurves(self, curve_list, color_list=["blue"], size_list = [0.5], custom_window=None, filename="Drifters"):
        ax = plt.axes(projection=ccrs.PlateCarree())
        ax.set_global()
        
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

        if self.storeFiles:
            plt.savefig(f"drifter_results/{int(self.points)}_{self.simpDELTA}_{self.freeDELTA}_{self.COMPLEXITY}_{filename}", dpi=300)
            plt.clf()
        else:
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
    # def getRossbyRadius(self, lat, lon):
    #     lon = (lon + 360) % 360
    #     lat = math.floor(lat) + 0.5
    #     lon = math.floor(lon) + 0.5
    #     return self.rossbyRadii.get((lat, lon), 1)
    
    def drawInputField(self, filename="InputField"):
        self.__drawVectorField(self.datacurves, "blue", filename)

    def drawResultField(self, relevant=False, filename="ResultField"):
        #interpolated_coords = [self.clustercurves[0]]  # Start with the first coordinate
        clustercurves = self.relevantClustercurves if relevant else self.clustercurves
        threshold = 100000
        new_curves = []
        for clustercurve in clustercurves:
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
        filename = filename + "Relevant" if relevant else filename+"All"
        self.__drawVectorField(new_curves, "red", filename)

    def __drawVectorField(self, curves, color, filename=""):
        
        lon, lat, u, v = [], [], [], []

        for path in curves:
            for i in range(len(path) - 1):
                x, y, z = path[i]
                lat1, lon1 = self.worldToLL(x, y, z)
                x, y, z = path[i+1]

                lat2, lon2 = self.worldToLL(x, y, z)
            
                delta_lon = lon2 - lon1
                delta_lat = lat2 - lat1
                mid_lon = (lon1 + lon2) / 2
                mid_lat = (lat1 + lat2) / 2
            
                lon.append(mid_lon)
                lat.append(mid_lat)
                u.append(delta_lon)
                v.append(delta_lat)

        lon = np.array(lon)
        lat = np.array(lat)
        u = np.array(u)
        v = np.array(v)

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
        #print(u_avg)

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
        if self.storeFiles:
            plt.savefig(f"drifter_results/{int(self.points)}_{self.simpDELTA}_{self.freeDELTA}_{self.COMPLEXITY}_{filename}", dpi=300)
            plt.clf()
        else:
            plt.show()

    
    
    def getOD(self):
        m = 0
        od = 0
        for curves, center in self.curves_per_cluster:
            print(".", end="")
            #center, curves 
            interpolated_curve = []
            
            for i in range(len(center) - 1):
                start_point = center[i]
                end_point = center[i + 1]
                
                segment_points = np.linspace(start_point, end_point, 5 + 2)
                
                interpolated_curve.append(segment_points)
            #print(center, interpolated_curve)
            if len(center) >= 2:
                center = np.vstack(interpolated_curve)
            for curve in curves:
                if len(curve) == 0 or len(center) == 0:
                    continue
                od += similaritymeasures.frechet_dist(center, curve)
                m += 1
        return od/m
    
    # def getSSE(self):
    #     sse = 0
    #     for cluster in self.curves_per_cluster:
    #         current_sse = 0
    #         print(".", end="")
    #         for c1 in cluster[0]:
    #             for c2 in cluster[0]:
    #                 current_sse += (similaritymeasures.frechet_dist(c1, c2))**2
    #         current_sse /= (2*len(cluster[0]) if len(cluster[0]) != 0 else 1)
    #         sse += current_sse
    #     return sse

    def getInitTime(self):
        return self.init_time
    
    def getExecutionTime(self):
        return self.execution_time
class KlClusterDriftersSolver(DriftersSolver):
    def __init__(self, drifterFiles, points=100000, simpDelta=20_000, freeDelta=200_000, complexity=1, iterations=0, storeFiles=False):
        super().__init__(drifterFiles, points)

        self.curves = kl.Curves()
        self.storeFiles = storeFiles
        self.iterations = iterations
        dc = []
        for curvedata in self.datacurves:
            lat, lon = self.worldToLL(curvedata[0][0], curvedata[0][1], curvedata[0][2])
            dc.append(curvedata)
        self.datacurves = dc

        for i, curvedata in enumerate(self.datacurves):
        #     #print(curvedata)
             curve = kl.Curve(curvedata, f"Curve_{i}")
        #     weights = curve.getWeights()
        #     if rossbyRadius:
        #         for j, p in enumerate(curve):
        #             lat, lon = self.worldToLL(p[0], p[1], p[2])
        #             weights[j] = (self.getRossbyRadius(lat, lon)/3)+10
        #             #if weights[j] != 500:
        #             #print(weights[j])
        #     curve.setWeights(weights)
             self.curves.add(curve)

        # self.DELTA = 50000
        # self.freeDELTA = 300000
        # self.simpDELTA = 5000
        self.freeDELTA = freeDelta
        self.simpDELTA = simpDelta
        self.COMPLEXITY = complexity
        self.ROUNDS = 1
        

        print(f"Inititalizing {len(self.curves)} curves")
        self.cc = kl.CurveClusterer()
        #if rossbyRadius:
        #    self.cc.initCurvesDiffDelta(self.curves, self.simpDELTA, self.freeDELTA)
        #else:
            #self.cc.initCurves(self.curves, self.DELTA)
        if self.iterations == 0:
            self.cc.initCurvesDiffDelta(self.curves, self.simpDELTA, self.freeDELTA)
        else:    
            self.init_time = (timeit.timeit(lambda: self.cc.initCurvesDiffDelta(self.curves, self.simpDELTA, self.freeDELTA), number=self.iterations) / self.iterations) if self.iterations != 0 else 0

            #self.cc.initCurvesDiffDelta(self.curves, self.simpDELTA, self.freeDELTA)

        
        self.simplifiedCurves = self.cc.getSimplifications()

    def solve(self, onlyRelevantClusters = False, withShow = False):
        #print(self.DELTA, self.curves)
        self.execution_time = (timeit.timeit(lambda: self.cc.greedyCover(self.COMPLEXITY, self.ROUNDS, withShow), number=self.iterations) / self.iterations) if self.iterations != 0 else 0

        self.clusters = self.cc.greedyCover(self.COMPLEXITY, self.ROUNDS, withShow)
        print(".")
        filterCount = 0
        self.curves = self.cc.getCurves()
        # for i in range(len(self.curves)):
        #     print(self.curves[i][0][0], self.simplifiedCurves[i][0][0], "????")

        # convert cluster centers to np array
        self.clustercurves = []
        self.relevantClustercurves = []

        #for i in range(len(self.curves)):
        #    print(self.curves[i][0][0], self.simplifiedCurves[i][0][0])



        self.curves_per_cluster = []
        x = 0
        for cluster in self.clusters:
            x+= 1
            print(cluster, len(self.clusters), x)
            curve = []


            center = cluster.center()
            start = int(center.start.value)
            end = int(center.end.value)
            curve = self.simplifiedCurves[center.curve]
            relevant=True
            if onlyRelevantClusters:
                if end-start <= 1 or len(cluster.values()) <= 1:
                    filterCount += 1
                    relevant = False
                    continue
            curvedata = []
            for i in range(start, end+1):
                if (i < 0 or i >= len(curve)):
                    raise Exception("Index out of bounds")
                curvedata.append(curve[i].values)

            all_curves = []
            for c in cluster:
                
                start = int(c.start.value)
                end = int(c.end.value)
                curve = self.curves[c.curve]
                curvedata_ = []
                start = int(self.cc.mapToBase(c.curve, c.start).value)
                end = int(self.cc.mapToBase(c.curve, c.end).value)
                for i in range(start, end+1):
                    if (i < 0 or i >= len(curve)):
                        raise Exception("Index out of bounds")
                    curvedata_.append(curve[i].values)
                all_curves.append(np.array(curvedata_))#self.simplifiedCurves[c.curve]))
            if relevant:
                self.relevantClustercurves.append(np.array(curvedata))
                self.relevantCurves_per_cluster.append((all_curves, np.array(curvedata)))
            self.clustercurves.append(np.array(curvedata))
            self.curves_per_cluster.append((all_curves, np.array(curvedata)))

        #self.curves = self.cc.getCurves()
        if onlyRelevantClusters:
            print(f"Filtered out {filterCount} center curves")
    
