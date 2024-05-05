import numpy as np
import matplotlib.pyplot as plt
import json
import klcluster as kl

GROUND_TRUTH_PATH = "../data/groundtruth86.json"


def getGroundTruth(instanceid):
    """get ground truth of cmu instance {instanceid>=1}"""
    labelToIndex = {}
    labelCount = 1

    def getLabelIndex(label):
        if not label in labelToIndex:
            nonlocal labelCount
            labelToIndex[label] = labelCount
            labelCount += 1
        return labelToIndex[label]

    gt = kl.GroundTruth()
    with open(GROUND_TRUTH_PATH) as f:
        jsonGt = json.loads(f.read())
        for i, row in enumerate(jsonGt["gt"][instanceid - 1]):
            gt.add(getLabelIndex(jsonGt["annotations"][instanceid - 1][i]), row[0])
            if row[0] != row[2]:
                gt.add(0, row[2])

        return gt


class CMUSolver:
    def __init__(self, filepath, groundTruth):
        print("Initializing " + filepath)
        trial = np.genfromtxt(filepath, delimiter=" ")
        print(trial.shape)

        self.curve = kl.Curve(trial)
        self.gt = groundTruth

        self.DELTA = 1.25
        self.COMPLEXITY = 10
        self.ROUNDS = 1
        self.cc = kl.CurveClusterer()

        curves = kl.Curves()
        curves.add(self.curve)
        groundTruths = kl.GroundTruths()
        groundTruths.add(self.gt)
        self.cc.initCurvesWithGT(curves, self.DELTA, groundTruths)

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

    def solve(self):
        return self.cc.greedyCover(self.COMPLEXITY, self.ROUNDS)

    # todo: label based on base curve
    def labelClustersBestMatch(self, clusters):
        labels = []
        gt = self.cc.getSimplifiedGTs()[0]
        # label all clusters
        for cluster in clusters:
            labelcount = {}
            # compute count for each label
            for matching in cluster:
                for i in range(int(matching.start.value), int(matching.end.value) + 1):
                    gtIndex = 0
                    while int(gt.paramAt(gtIndex).value) < i and gtIndex < len(gt):
                        gtIndex += 1
                    label = gt.labelAt(gtIndex)
                    if not label in labelcount:
                        labelcount[label] = 0
                    labelcount[label] += 1
            # pick maximum
            if len(labelcount) > 0:
                bestLabel = max(labelcount, key=labelcount.get)
            else:
                bestLabel = -1
            labels.append(bestLabel)

        return labels

    def getBaseSementation(self, clusters):
        labels = self.labelClustersBestMatch(clusters)

        size = self.curve.complexity
        segmentation = [-1] * size
        for clusterIdx, cluster in enumerate(clusters):
            if labels[clusterIdx] == -1:
                continue
            for matching in cluster:
                start = int(self.cc.mapToBase(0, matching.start).value)
                end = int(self.cc.mapToBase(0, matching.end).value)

                # print(f"{start} bis {end} : {labels[clusterIdx]}")

                for i in range(start, end+1):
                    if (i < 0 or i >= size):
                        raise Exception("Index out of bounds")

                    if segmentation[i] == -1 or segmentation[i] == labels[clusterIdx]:
                        segmentation[i] = labels[clusterIdx]
                    else:
                        segmentation[i] = 0

        return segmentation

    def getSegments(self, segmentation):
        segments = []
        segmentStart = 0

        for i in range(1, len(segmentation)):
            if segmentation[i] != segmentation[i - 1]:
                segment = {"size": i - segmentStart, "label": segmentation[i - 1]}
                segments.append(segment)
                segmentStart = i

        # append last segment
        if len(segmentation) > 0:
            segment = {"size": len(segmentation) - segmentStart, "label": segmentation[i - 1]}
            segments.append(segment)

        return segments

    def getGTSegments(self):
        segments = []
        segmentStart = 0
        for i in range(0, len(self.gt)):
            segment = {"size": self.gt[i][1] - segmentStart, "label": self.gt[i][0]}
            segments.append(segment)
            segmentStart = self.gt[i][1] + 1
        return segments
    
    def plotSegmentation(self, segmentation):
        segments = self.getSegments(segmentation)
    
        plt.rcParams['toolbar'] = 'None' # Remove tool bar
        fig, ax = plt.subplots(tight_layout=True)
        fig.canvas.manager.set_window_title("Curve segmentation")
        fig.set_size_inches(10, 1.5)
        self.__plotSegmentsToSubplot(ax, segments)
        plt.show()

    def plotSegmentationAndGT(self, segmentation):
        segments = self.getSegments(segmentation)
        gtSegments = self.getGTSegments()
    
        plt.rcParams['toolbar'] = 'None' # Remove tool bar
        fig, axs = plt.subplots(2, 1, tight_layout=True)
        fig.canvas.manager.set_window_title("Curve segmentation compared to GT")
        fig.set_size_inches(10, 3)
        axs[0].set_title("GT")
        self.__plotSegmentsToSubplot(axs[0], gtSegments)
        axs[1].set_title("Clustering")
        self.__plotSegmentsToSubplot(axs[1], segments)
        plt.show()

    def __plotSegmentsToSubplot(self, ax, segments):
        x = 0
        for s in segments:
            if s["label"] < len(self.colors):
                ax.barh(0, s["size"], left=x, color=self.colors[s["label"]])
            else:
                ax.barh(0, s["size"], left=x, label=s["label"])

            x += s["size"]

        xax = ax.axes.get_yaxis()
        xax.set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_visible(False)
        
        

solver = CMUSolver("../data/86_1.txt", getGroundTruth(1))
result = solver.solve()

segmentation = solver.getBaseSementation(result)

solver.plotSegmentationAndGT(segmentation)
