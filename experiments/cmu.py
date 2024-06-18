import numpy as np
import matplotlib.pyplot as plt
import json
import klcluster as kl

GROUND_TRUTH_PATH = "../data_cmu/groundtruth86.json"


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
    
    def labelClustersBestMatch(self, clusters):
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

    def getBaseSementation(self, clusters):
        labels = self.labelClustersBestMatch(clusters)

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
                segment = {"size": i - segmentStart, "label": segmentation[i - 1]}
                segments.append(segment)
                segmentStart = i

        # append last segment
        if len(segmentation) > 0:
            segment = {"size": len(segmentation) - segmentStart, "label": segmentation[i - 1]}
            segments.append(segment)

        return segmentation, segments

    def getGTSegments(self):
        segments = []
        segmentStart = 0
        for i in range(0, len(self.gt)):
            segment = {"size": self.gt[i][1] - segmentStart, "label": self.gt[i][0]}
            segments.append(segment)
            segmentStart = self.gt[i][1] + 1
        print(segments)
        return segments
    
    def plotSegmentation(self, segments):    
        plt.rcParams['toolbar'] = 'None' # Remove tool bar
        fig, ax = plt.subplots(tight_layout=True)
        fig.canvas.manager.set_window_title("Curve segmentation")
        fig.set_size_inches(10, 1.5)
        self.__plotSegmentsToSubplot(ax, segments)
        plt.show()

    def plotSegmentationAndGT(self, segments):
        gtSegments = self.getGTSegments()
        print(segments)
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
                ax.barh(0, s["size"], left=x, edgecolor="black", linewidth=1.5, color=self.colors[s["label"]])
            else:
                ax.barh(0, s["size"], left=x, edgecolor="black", linewidth=1.5, label=s["label"])

            x += s["size"]

        xax = ax.axes.get_yaxis()
        xax.set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_visible(False)
    
    def calculateAccurcacy(self, segmentation):
        M = self.getConfusionMatrix(segmentation)
        # main diagonal divided my all entries of the confusion matrix
        return sum([M[i][i] for i in range(len(M))]) / sum([M[i][j] for i in range(len(M)) for j in range(len(M))])

    def calculateClassPrecision(self, segmentation, classLabel):
        M = self.getConfusionMatrix(segmentation)
        # main diagonal divided my all entries of the confusion matrix
        positiveClassified = sum([M[classLabel][i] for i in range(len(M))])
        return M[classLabel][classLabel] / positiveClassified

    def calculateMacroPrecision(self, segmentation):
        M = self.getConfusionMatrix(segmentation)
        return sum([self.calculateClassPrecision(segmentation, label) for label in range(len(M))]) / len(M)


    def calculateClassRecall(self, segmentation, classLabel):
        M = self.getConfusionMatrix(segmentation)
        # main diagonal divided my all entries of the confusion matrix
        GTPositives = sum([M[i][classLabel] for i in range(len(M))])
        return M[classLabel][classLabel] / GTPositives

    def calculateMacroRecall(self, segmentation):
        M = self.getConfusionMatrix(segmentation)
        return sum([self.calculateClassRecall(segmentation, label) for label in range(len(M))]) / len(M)

    def getConfusionMatrix(self, segmentation):
        gtSegments = self.getGTSegments()
        num_labels = max([gtSegment['label'] for gtSegment in gtSegments])
        M = [[0 for i in range(5)] for j in range(5)]
        idx = 0
        for gtSegment in gtSegments:
            for i in range(idx, idx+gtSegment['size']):
                M[segmentation[i]][gtSegment['label']] +=1
            idx += gtSegment['size']
        return M

solver = CMUSolver("../data_cmu/86_1.txt", getGroundTruth(1))
print(solver.gt)
result = solver.solve()

segmentation, segments = solver.getBaseSementation(result)
print(segmentation)
solver.plotSegmentationAndGT(segments)
print(solver.calculateAccurcacy(segmentation), solver.calculateMacroPrecision(segmentation), solver.calculateMacroRecall(segmentation))