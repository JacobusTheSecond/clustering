from abc import ABC, abstractmethod
import numpy as np
import matplotlib.pyplot as plt
import csv
import os
import json
from dataclasses import dataclass
from collections import Counter
import timeit

import matlab.engine
import klcluster as kl
import csv
import math



class IMUSolver(ABC):
    def __init__(self):
        
        self.DATA_PATH = os.path.join(os.path.dirname(__file__), "../xMFs_data/")
        
        self.instances = [2, 3, 4, 5, 6, 8, 9]
        self.activities = ["flat", "stair", "stop"]
        self.curveData = []
        self.gtData = []
        self.fullgtData = []

        for inst in self.instances:
            for a in self.activities:
                path = self.DATA_PATH + f"xMF0{inst}/{a}"
                for subdir in os.listdir(path):
                        
                    with open(path + f"/{subdir}/norm/jointDataNorm.csv", 'r', newline='') as csvfile:
                        csvreader = csv.reader(csvfile)
                        
                        selected_lines = []
                        
                        # Loop through the CSV file
                        for i, row in enumerate(csvreader):
                            # Break if we've processed 100 lines
                            if len(selected_lines) >= 100:
                                break
                            # Only process every 5th line
                            if (i + 1) % 5 == 0:
                                row = row[3:]
                                for i in range(len(row)):
                                    row[i] = float(row[i])
                                    # normalization for data has been formatted wrong
                                    if row[i] > 1:
                                        row[i] /= 10 ** math.floor(math.log10(row[i]))
                                selected_lines.append(row)
                    
                    curve = np.array(selected_lines)
                    self.curveData.append(curve)
                    self.gtData.append(a)
                    self.fullgtData.append([inst, a, subdir])
        

    @abstractmethod
    def solve():
        pass


class KlClusterIMUSolver(IMUSolver):
    def __init__(self, SIMP_DELTA = 1.25, FREE_DELTA = 1.25, COMPLEXITY = 10):
        super().__init__()

        map_label = {"flat" : 0,
                     "stair" : 1,
                     "stop" : 2} 

        self.curves = kl.Curves()
        self.SIMP_DELTA = SIMP_DELTA
        
        self.groundTruths = kl.GroundTruths()

        for c, label in zip(self.curveData, self.gtData):
            curve = kl.Curve(c, f"Curve")
            gt = kl.GroundTruth()
            for i in range(101):
                gt.add(map_label[label], i)
                
            print(len(gt), len(curve))
            self.curves.add(curve)
            self.groundTruths.add(gt)

        print(len(self.groundTruths), len(self.curves))
        self.cc = kl.CurveClusterer()
        self.cc.initCurvesWithGTDiffDelta(self.curves, 0.8, 1.7, self.groundTruths)
        # self.cc.initCurvesWithGTDiffDelta(self.curves, 0.7, 2.3, self.groundTruths)


    def solve(self):
        self.clusters = self.cc.greedyCover(3, 1, False)
        print(len(self.clusters))

        m = 0
                
        for cluster in self.clusters:
            labelcount = {}
            totalcount = 0
            print(cluster)
            # compute count for each label
            for matching in cluster:
                #print(matching.curve)
                #print(matching.values()[0], m)
                matchingstart = int(self.cc.mapToBase(0, matching.start).value)
                matchingend = int(self.cc.mapToBase(0, matching.end).value)
                if matchingend - matchingstart >= 10:
                    #print(matchingstart, matchingend)
                    print(matchingstart, matchingend, matching.curve)
                    print(matching.values(), matching.values().shape)
                    #print(self.gtData[int(matching.values()[0])], self.fullgtData[int(matching.values()[0])], matchingstart, matchingend)

                # matchingstart = int(self.cc.mapToBase(0, matching.start).value)
                # matchingend = int(self.cc.mapToBase(0, matching.end).value)
                # print(matchingstart, matchingend)