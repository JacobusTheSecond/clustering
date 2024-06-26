from abc import ABC, abstractmethod
import numpy as np
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import csv

import klcluster as kl

class DriftersSolver(ABC):
    def __init__(self, drifterFiles):
        self.clustercurves = None
        self.datacurves = []
        for filepath in drifterFiles:
            with open(filepath, "r") as f:
                curveData = np.array(list(csv.reader(f, delimiter=" "))).astype(float)
                self.datacurves.append(curveData)

    @abstractmethod
    def solve():
        pass
    
    def plotInput(self):
        self.__plotCurves(self.datacurves)

    def plotResults(self):
        if (not self.clustercurves):
            raise Exception("No result found, consider calling solve() before")
        self.__plotCurves(self.clustercurves, "red")

    def plotInputAndResult(self):
        self.__plotCurves(self.datacurves, "blue", self.clustercurves, "red")

    def __plotCurves(self, curves, color="blue", curves2=None, color2=None):
        ax = plt.axes(projection=ccrs.PlateCarree())
        ax.set_global()
        #ax.coastlines()
        ax.stock_img()
        plt.tight_layout()
        ax.get_figure().canvas.manager.set_window_title("Drifters")

        def add_arrow(line, size=0.001, color=None):
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
            plt.arrow(x, y, dx, dy, facecolor=color, edgecolor=color,
            head_width=1, head_length=1)

            # line.axes.annotate('',
            #     [x + dx, y + dy],
            #     [x, y],
            #     arrowprops=dict(width=0, headwidth=3, headlength=3, color=color),
            #     annotation_clip=True
            # )

        for k in range(0, 2):
            if curves == None:
                break
            for i, curve in enumerate(curves):
                if curve.size < 3:
                    continue
                lat, lon = self.worldToLL(curve[:,0], curve[:,1], curve[:,2])

                line = plt.plot(lon, lat,
                        color=color, linewidth=0.7,
                        transform=ccrs.Geodetic())[0]
                
                add_arrow(line)

                if i % 100 == 0:
                    print(f"{i}/{len(curves)}")
                i += 1
                # if i > 600:
                #     break

            curves = curves2 # switch to second curve set
            color = color2

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
    

class KlClusterDriftersSolver(DriftersSolver):
    def __init__(self, drifterFiles):
        super().__init__(drifterFiles)

        self.curves = kl.Curves()
        for i, curvedata in enumerate(self.datacurves):
            self.curves.add(kl.Curve(curvedata, f"Curve_{i}"))

        self.DELTA = 15000
        self.COMPLEXITY = 10
        self.ROUNDS = 1

        self.cc = kl.CurveClusterer()
        self.cc.initCurves(self.curves, self.DELTA)
        self.simplifiedCurves = self.cc.getSimplifications()

    def solve(self):
        self.clusters = self.cc.greedyCover(self.COMPLEXITY, self.ROUNDS)

        # convert cluster centers to np array
        self.clustercurves = []
        for cluster in self.clusters:
            curve = []

            center = cluster.center()
            start = int(center.start.value)
            end = int(center.end.value)
            curve = self.simplifiedCurves[center.curve]

            curvedata = []
            for i in range(start, end+1):
                if (i < 0 or i >= len(curve)):
                    raise Exception("Index out of bounds")
                # print(curve[i])
                curvedata.append(curve[i].values)

            self.clustercurves.append(np.array(curvedata))
