from abc import ABC, abstractmethod
import numpy as np
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import csv
import os
import json
from dataclasses import dataclass

import klcluster as kl

@dataclass
class Segment:
    size: str
    label: float

class CMUSolver(ABC):
    def __init__(self, tag):
        if tag < 1 or tag > 14:
            raise Exception("CMU trial tag must be between 1 and 14")
        self.DATA_PATH = os.path.join(os.path.dirname(__file__), "../data_cmu/")
        self.tag = tag
        self.gtSegments = self.__getGTSegments(tag)

        self.segmentation = None
        self.segments = None

        self.colors = [
            (0, 0, 0),
            (1, 1, 0),
            (1, 0.75, 0.25),
            (1, 0.5, 0.5),
            (1, 0.25, 0.75),
            (1, 0, 1),
            (0.75, 0.25, 1),
            (0.5, 0.5, 1),
            (0.25, 0.75, 1),
            (0, 1, 1),
            (0.25, 1, 0.75),
            (0.5, 1, 0.5),
            (0.75, 1, 0.25),
        ]

    @abstractmethod
    def solve():
        pass

    def plotSegments(self, segments):    
        plt.rcParams['toolbar'] = 'None' # Remove tool bar
        fig, ax = plt.subplots(tight_layout=True)
        fig.canvas.manager.set_window_title("Curve segmentation")
        fig.set_size_inches(10, 1.5)
        self.__plotSegmentsToSubplot(ax, segments)
        plt.show()

    def plotSegmentsAndGT(self, segments):    
        plt.rcParams['toolbar'] = 'None' # Remove tool bar
        fig, axs = plt.subplots(2, 1, tight_layout=True)
        fig.canvas.manager.set_window_title("Curve segmentation compared to GT")
        fig.set_size_inches(10, 3)
        axs[0].set_title("GT")
        self.__plotSegmentsToSubplot(axs[0], self.gtSegments)
        axs[1].set_title("Clustering")
        self.__plotSegmentsToSubplot(axs[1], segments)
        plt.show()

    def __plotSegmentsToSubplot(self, ax, segments):
        x = 0
        for s in segments:
            if s.label < len(self.colors):
                ax.barh(0, s.size, left=x, edgecolor="black", linewidth=1.5, color=self.colors[s.label])
            else:
                ax.barh(0, s.size, left=x, edgecolor="black", linewidth=1.5, label=s.label)

            x += s.size

        xax = ax.axes.get_yaxis()
        xax.set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_visible(False)
    
    def __getGTSegments(self, tag):
        """get ground truth of cmu instance {instanceid>=1}"""
        labelToIndex = {}
        labelCount = 1

        def getLabelIndex(label):
            if not label in labelToIndex:
                nonlocal labelCount
                labelToIndex[label] = labelCount
                labelCount += 1
            return labelToIndex[label]

        segments = []
        segmentStart = 0
        with open(os.path.join(self.DATA_PATH, "groundtruth86.json")) as f:
            jsonGt = json.loads(f.read())
            for i, row in enumerate(jsonGt["gt"][tag - 1]):
                # gt.add(getLabelIndex(jsonGt["annotations"][tag - 1][i]), row[0])
                segments.append(Segment(size=row[0] - segmentStart, label=getLabelIndex(jsonGt["annotations"][tag - 1][i])))
                segmentStart = row[0]
                if row[0] != row[2]:
                    segments.append(Segment(size=row[2] - segmentStart, label=0))
                    segmentStart = row[2]
                    # gt.add(0, row[2])

        return segments

class KlClusterCMUSolver(CMUSolver):
    def __init__(self, tag):
        super().__init__(tag)

        trial = np.genfromtxt(os.path.join(self.DATA_PATH, f"86_{tag}.txt"), delimiter=" ")

        self.curve = kl.Curve(trial, f"trial_{tag}")
        curves = kl.Curves()
        curves.add(self.curve)

        # construct gt
        self.gt = kl.GroundTruth()
        index = 0
        for segment in self.gtSegments:
            index += segment.size
            self.gt.add(segment.label, index)
        groundTruths = kl.GroundTruths()
        groundTruths.add(self.gt)
        
        self.DELTA = 1.25
        self.COMPLEXITY = 10
        self.ROUNDS = 1
        self.cc = kl.CurveClusterer()
        self.cc.initCurvesWithGT(curves, self.DELTA, groundTruths)

    def solve(self):
        clusters = self.cc.greedyCover(self.COMPLEXITY, self.ROUNDS)
        self.segments, self.segmentation = self.__getBaseSementation(clusters)
        return self.segments, self.segmentation

    def __getBaseSementation(self, clusters):
        labels = self.__labelClustersBestMatch(clusters)

        # build segmentation
        size = self.curve.complexity
        segmentation = [-1] * size
        is_matching_startend = [False] * size # store whether curve point is start or end point of a cluster matching
        for clusterIdx, cluster in enumerate(clusters):
            if labels[clusterIdx] == -1:
                continue
            for matching in cluster:
                start = int(self.cc.mapToBase(0, matching.start).value)
                end = int(self.cc.mapToBase(0, matching.end).value)

                is_matching_startend[start] = True
                is_matching_startend[end] = True

                # print(f"{start} bis {end} : {labels[clusterIdx]}")

                for i in range(start, end+1):
                    if (i < 0 or i >= size):
                        raise Exception("Index out of bounds")
                    if segmentation[i] == -1 or segmentation[i] == labels[clusterIdx]:
                        segmentation[i] = labels[clusterIdx]
                    else:
                        segmentation[i] = 0

        # turn all unknown (-1) to black (0)
        for i in range(0, len(segmentation)):
            if segmentation[i] == -1:
                segmentation[i] = 0

        # turn segmentation to segments
        segments = []
        segmentStart = 0

        for i in range(1, len(segmentation)):
            if segmentation[i] != segmentation[i - 1] or is_matching_startend[i]:
                segments.append(Segment(size=i - segmentStart, label=segmentation[i - 1]))
                segmentStart = i

        # append last segment
        if len(segmentation) > 0:
            segments.append(Segment(size=len(segmentation) - segmentStart, label=segmentation[i - 1]))

        return segments, segmentation
    
    def __labelClustersBestMatch(self, clusters):
        labels = []
        # label all clusters
        for cluster in clusters:
            labelcount = {}
            totalcount = 0
            # compute count for each label
            for matching in cluster:
                matchingstart = int(self.cc.mapToBase(0, matching.start).value)
                matchingend = int(self.cc.mapToBase(0, matching.end).value)
                for i in range(matchingstart, matchingend + 1):
                    gtIndex = 0
                    while gtIndex < len(self.gt) and self.gt[gtIndex][1] < i:
                        gtIndex += 1
                    if gtIndex >= len(self.gt):
                        raise Exception("Matching index out of bounds of GT")
                    label = self.gt[gtIndex][0]
                    if not label in labelcount:
                        labelcount[label] = 0
                    labelcount[label] += 1
                    totalcount += 1
            # pick maximum
            if len(labelcount) > 0:
                bestLabel = max(labelcount, key=labelcount.get)
                # if label has no majority, discard cluster
                if labelcount[bestLabel] * 2 <= totalcount:
                    bestLabel = -1
            else:
                bestLabel = -1
            labels.append(bestLabel)

        return labels
