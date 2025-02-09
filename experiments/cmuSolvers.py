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
        self.gtSegments, self.labelMapping = self.getGTSegments(tag)

        self.frames = sum([int(segment.size) for segment in self.gtSegments])
        self.segmentation = None
        self.segments = None
        self.execution_time = None

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

    def solveAndTime(self, timing_rounds=0, *args, **kwargs):
        self.execution_time = (timeit.timeit(lambda: self.solve(*args, **kwargs), number=timing_rounds) / timing_rounds) if timing_rounds != 0 else 0
        return self.solve()

    def plotSegments(self, segmentsList, methodNames=None):
        plt.rcParams['toolbar'] = 'None' # Remove tool bar
        fig, axs = plt.subplots(len(segmentsList), 1, tight_layout=True, squeeze=False)
        fig.canvas.manager.set_window_title("Curve segmentation")
        fig.set_size_inches(10, len(segmentsList) * 1.5)
        for i in range(0, len(segmentsList)):
            if methodNames and i < len(methodNames):
                axs[i][0].set_title(methodNames[i])
            self.plotSegmentsToSubplot(axs[i][0], segmentsList[i])
        plt.show()

    def plotSegmentsAndGT(self, segmentsList, methodNames=None):
        plt.rcParams['toolbar'] = 'None' # Remove tool bar
        fig, axs = plt.subplots(len(segmentsList)+1, 1, tight_layout=True, squeeze=False)
        fig.canvas.manager.set_window_title("Curve segmentation compared to GT")
        fig.set_size_inches(10, (len(segmentsList)+1) * 1.5)
        axs[0][0].set_title("GT")
        self.plotSegmentsToSubplot(axs[0][0], self.gtSegments)
        for i in range(0, len(segmentsList)):
            if methodNames and i < len(methodNames):
                axs[i+1][0].set_title(methodNames[i])
            self.plotSegmentsToSubplot(axs[i+1][0], segmentsList[i])
        plt.savefig("overview_acc.png")
        plt.show()
    
    def plotSegmentsToSubplot(self, ax, segments):
        x = 0
        for s in segments:
            if s.label < len(self.colors):
                ax.barh(0, s.size, left=x, edgecolor="grey", linewidth=0.5, color=self.colors[s.label])
            else:
                ax.barh(0, s.size, left=x, edgecolor="grey", linewidth=1.5, label=s.label)

            x += s.size

        xax = ax.axes.get_yaxis()
        xax.set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_visible(False)
    
    def getGTSegments(self, tag):
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

        return segments, labelToIndex
    
    def calculateAccuracy(self, segmentation, tag):
        M = self.getConfusionMatrix(segmentation, tag)
        # main diagonal divided my all entries of the confusion matrix

        return sum([M[i][i] for i in range(len(M))]) / sum([M[i][j] for i in range(len(M)) for j in range(len(M))])

    def calculateAccuracyTrueLabels(self, segmentation, tag):
        M = self.getConfusionMatrix(segmentation, tag)
        return sum([M[i][i] for i in range(1, len(M))]) / sum([M[i][j] for i in range(0, len(M)) for j in range(1, len(M))])


    def calculateClassPrecision(self, segmentation, classLabel, tag):
        M = self.getConfusionMatrix(segmentation, tag)
        positiveClassified = sum([M[classLabel][i] for i in range(len(M))])
        return M[classLabel][classLabel] / positiveClassified if positiveClassified != 0 else 0

    def calculateMacroPrecision(self, segmentation, tag):
        M = self.getConfusionMatrix(segmentation, tag)
        return sum([self.calculateClassPrecision(segmentation, label, tag) for label in range(len(M))]) / len(M)


    def calculateClassRecall(self, segmentation, classLabel, tag):
        M = self.getConfusionMatrix(segmentation, tag)
        GTPositives = sum([M[i][classLabel] for i in range(len(M))])
        return M[classLabel][classLabel] / GTPositives if GTPositives != 0 else 1

    def calculateMacroRecall(self, segmentation, tag):
        M = self.getConfusionMatrix(segmentation, tag)
        return sum([self.calculateClassRecall(segmentation, label, tag) for label in range(len(M))]) / len(M)

    def calculateMacroF1(self, segmentation, tag):
        prec = self.calculateMacroPrecision(segmentation, tag)
        recall = self.calculateMacroRecall(segmentation, tag)
        return (2*prec*recall) / (prec + recall)

    def getConfusionMatrix(self, segmentation, tag):
        gtSegments = self.getGTSegments(tag)
        num_labels = max([gtSegment.label for gtSegment in gtSegments[0]] + [segment.label for segment in segmentation])+1
        M = [[0 for i in range(num_labels)] for j in range(num_labels)]
        frameListGT = []
        frameList = []
        for gtSegment in gtSegments[0]:
            frameListGT.extend([gtSegment.label]*gtSegment.size)
        for segment in segmentation:
            frameList.extend([segment.label]*segment.size)
        for frame, frameGT in zip(frameList, frameListGT):
            M[frame][frameGT] +=1
        return M
    
    def getExecutionTime(self):
        return self.execution_time

class KlClusterCMUSolver(CMUSolver):
    def __init__(self, tag, SIMP_DELTA = 1.25, FREE_DELTA = 1.25, COMPLEXITY = 10):
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
        
        self.SIMP_DELTA = SIMP_DELTA
        self.FREE_DELTA = FREE_DELTA
        self.COMPLEXITY = COMPLEXITY
        self.ROUNDS = 1
        self.cc = kl.CurveClusterer()
        self.cc.initCurvesWithGTDiffDelta(curves, self.SIMP_DELTA, self.FREE_DELTA, groundTruths)

    def solve(self, mergeOverlappingClusters = False):
        #self.execution_time = (timeit.timeit(lambda: self.cc.greedyCover(self.COMPLEXITY, self.ROUNDS, False), number=timing_rounds) / timing_rounds) if timing_rounds != 0 else 0
        #print("EXECUTIN TIME", self.execution_time, timing_rounds)
        self.clusters = self.cc.greedyCover(self.COMPLEXITY, self.ROUNDS, False)
        if mergeOverlappingClusters:
            print(f"Merged {self.cc.mergeOverlappingClusters(self.clusters, 0.5)} clusters.")
        self.segments, self.segmentation = self.__getBaseSementation(self.clusters)
        return self.segments

    def __getBaseSementation(self, clusters):
        labels = self.__labelClustersBestMatch(clusters)
        self.clusterLabels = labels

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
    
    def __getBaseSementationMajority(self, clusters):
        labels = self.__labelClustersBestMatch(clusters)
        self.clusterLabels = labels

        # build segmentation
        size = self.curve.complexity
        segmentation = [[] for _ in range(size)]
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
                    segmentation[i].append(labels[clusterIdx])


        # compute majority label for each time step
        for i in range(0, len(segmentation)):
            if len(segmentation[i]) > 0:
                c = Counter(segmentation[i])
                most_common = c.most_common(1)[0]
                if most_common[1] > len(segmentation[i]) / 2: # check if label has majority
                    segmentation[i] = most_common[0]
                else:
                    segmentation[i] = 0
            else:
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

        # clusterGroups = []
        # groupIndexOfCluster = [i for i in range(len(clusters))]
        
        # for i in range(0, len(clusters)):
        #     for j in range(i + 1, len(clusters)):
        #         # calculate similarity
        #         overlap = 0

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
    
    def plotSegmentsAndMatching(self):
        plt.rcParams['toolbar'] = 'None' # Remove tool bar
        fig, axs = plt.subplots(3, 1, tight_layout=True, squeeze=False)
        fig.canvas.manager.set_window_title("Curve segmentation compared to GT")
        fig.set_size_inches(10, 3 * 1.5)
        axs[1][0].set_title("KlCluster")
        self.plotSegmentsToSubplot(axs[1][0], self.segments)
        axs[2][0].set_title("GT")
        self.plotSegmentsToSubplot(axs[2][0], self.gtSegments)
        
        # plot matchings of each cluster
        ax = axs[0][0]
        ax.set_title("Cluster matchings")
        for clusterIdx, cluster in enumerate(self.clusters):
            labelOrig = self.clusterLabels[clusterIdx]
            label = labelOrig
            if label == -1:
                label = 0
            for matching in cluster:
                matchingstart = int(self.cc.mapToBase(0, matching.start).value)
                matchingend = int(self.cc.mapToBase(0, matching.end).value)
                size = matchingend - matchingstart

                if labelOrig == -1:
                    ax.barh(clusterIdx, size, left=matchingstart, edgecolor="black", linewidth=0.7, color="gray")
                elif label < len(self.colors):
                    ax.barh(clusterIdx, size, left=matchingstart, edgecolor="black", linewidth=0.7, color=self.colors[label])
                else:
                    ax.barh(clusterIdx, size, left=matchingstart, edgecolor="black", linewidth=0.7, label=label)


        xax = ax.axes.get_xaxis()
        xax.set_visible(False)
        xax = ax.axes.get_yaxis()
        xax.set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.spines['bottom'].set_visible(False)

        plt.savefig(f"accuracy15.png")
        plt.show()


class AcaCMUSolver(CMUSolver):
    def __init__(self, tag, method):
        """ method: 'gmm' | 'aca' | 'haca' """
        super().__init__(tag)
        self.method = method

    def solve(self):
        print("start matlab")
        eng = matlab.engine.start_matlab()  # connect_matlab()
        print("finished")

        eng.cd(os.path.join(os.path.dirname(__file__), "../aca/aca"), nargout=0)
        #eng.make(nargout=0)
        eng.addPath(nargout=0) # add aca paths like src or lib
        eng.addpath("..") # add runAca, runGmm, ... to path

        match self.method:
            case "gmm":
                gt, seg, labelNames = eng.runGmm(self.tag, nargout=3)
                #self.execution_time = (timeit.timeit(lambda: eng.runGmm(self.tag, nargout=3), number=timing_rounds) / timing_rounds) if timing_rounds != 0 else 0

            case "aca":
                gt, seg, labelNames = eng.runAca(self.tag, nargout=3)
                #self.execution_time = (timeit.timeit(lambda: eng.runAca(self.tag, nargout=3), number=timing_rounds) / timing_rounds) if timing_rounds != 0 else 0

            case "haca":
                gt, seg, labelNames = eng.runHaca(self.tag, nargout=3)
                #self.execution_time = (timeit.timeit(lambda: eng.runHaca(self.tag, nargout=3), number=timing_rounds) / timing_rounds) if timing_rounds != 0 else 0

            case _:
                raise Exception(f"Unknown method '{self.method}'")
            

        eng.quit() # stop matlab.engine

        segFrames = np.array(seg["s"]).astype(int)[0] # convert matlab double array to numpy int array
        segLabels = np.array(seg["G"]).astype(int)

        # build segments
        segments = []
        segmentStart = segFrames[0] - 1
        for i in range(1, len(segFrames)):
            frame = segFrames[i] - 1
            labelIndices = np.argwhere(segLabels[:, i-1])
            if (len(labelIndices) > 1):
                raise Exception(f"More than one label for matlab segment {i}")
            label = self.labelMapping[labelNames[labelIndices.item()]] # map matlab label to our labels
            segments.append(Segment(size=frame - segmentStart, label=label))
            segmentStart = frame
        self.clusters = segments
        return segments


class TmmCMUSolver(CMUSolver):
    def __init__(self, tag):
        super().__init__(tag)
        
    def solve(self):
        print("start matlab")
        
        eng = matlab.engine.start_matlab()  # connect_matlab()
        print("finished")
        
        eng.cd(os.path.join(os.path.dirname(__file__), "..//MotionSegmentationCode/MatlabCodeTMM/"), nargout=0)
        # eng.make(nargout=0)
        eng.addPath("MotionSegmentation", nargout=0) # add tmm paths
        
        eng.cd("MotionSegmentation")
        #self.execution_time = (timeit.timeit(lambda: eng.call_segmentation(self.tag, nargout=3), number=timing_rounds) / timing_rounds) if timing_rounds != 0 else 0
        comps, sframes, eframes = eng.call_segmentation(self.tag, nargout=3)
        eng.quit() # stop matlab.engine

        comps = [int(i[0]) for i in comps]
        self.reducedFrames = int(eframes[-1][0])

        sframes = [int(self.frames/self.reducedFrames*(int(i[0])-1))+1 for i in sframes]
        eframes = [int(self.frames/self.reducedFrames*(int(i[0])-1))+1 for i in eframes]
        #eframes = [int(i[0]) for i in eframes]
        
        #segments_ = [[, int(self.frames/self.reducedFrames*(end-1))+1, label] for start, end, label in segments_]
        

        self.clusters = [[] for i in range(max(comps))]
        ## stop matlab.engine
        for i in range(len(comps)):
            self.clusters[comps[i]-1].append((sframes[i], eframes[i]))

        self.gt = self.getGTSegments(self.tag)
        frameListGT = []
        for gtSegment in self.gt[0]:
            frameListGT.extend([gtSegment.label]*gtSegment.size)
        segments_ = []
        labels = []
        for clusterIdx, cluster in enumerate(self.clusters):
            labelcount = {}
            totalcount = 0
            for subsegment in cluster:
                start = subsegment[0]
                end = subsegment[1]
                
                for i in range(start-1, end - 1):
                   
                    label = frameListGT[i]
                    if not label in labelcount:
                        labelcount[label] = 0
                    labelcount[label] += 1
                    totalcount += 1
            # pick maximum
            if len(labelcount) > 0:
                bestLabel = max(labelcount, key=labelcount.get)
                # if label has no majority, discard cluster
                if labelcount[bestLabel] * 2 <= totalcount:
                    bestLabel = 0
            else:
                bestLabel = 0
            labels.append(bestLabel)
            for subsegment in cluster:
                start = subsegment[0]
                end = subsegment[1]
                segments_.append([start, end, bestLabel])
        segments_.sort()
        self.reducedFrames = segments_[-1][1]
       # segments_ = [[int(self.frames/self.reducedFrames*(start-1))+1, int(self.frames/self.reducedFrames*(end-1))+1, label] for start, end, label in segments_]
        segments = [Segment(size=s[1]-s[0], label=s[2]) for s in segments_]
        print(sum([int(segment.size) for segment in segments]))

        print("Anzahl Labels Tmm", len(labels), len(self.clusters), len(segments))
        return segments
