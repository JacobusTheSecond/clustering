from cmuSolvers import KlClusterCMUSolver
from cmuSolvers import AcaCMUSolver
from cmuSolvers import TmmCMUSolver
import pandas as pd
import os
import timeit
import json

file_path = "result.json"

result = {}

ITERATIONS = 3

columns = ["TAG", "Method", "nClusters", "Time", "Accuracy", "True Accuracy", "Macro Precision", "Macro Recall", "Macro F1", "Num iterations", "Complexity", "Simp Delta", "Free Delta"]

csv_file = 'outputCMU.csv'

if os.path.exists(csv_file) and os.path.getsize(csv_file) > 0:
    df = pd.read_csv(csv_file)

    if df.empty:
        df = pd.DataFrame(columns=columns)
else:
    df = pd.DataFrame(columns=columns)


configurations = [
    [[TAG for TAG in range(1, 2)], [["KlCluster", 19, 0.8, 1.1]], 3],
                [[TAG for TAG in range(1, 2)], [["Aca"], ["Haca"], ["Gmm"]], 3],
                [[TAG for TAG in range(1, 2)], [["Tmm"]], 3],

                [[TAG for TAG in range(1, 15)], [
                    ["KlClusterS", 5, 0.2, 1],
                    ["KlClusterS", 10, 0.2, 1],
                    ["KlClusterS", 15, 0.2, 1],
                    ["KlClusterS", 5, 0.3, 1],
                    ["KlClusterS", 10, 0.3, 1],
                    ["KlClusterS", 15, 0.3, 1],
                    ["KlClusterS", 5, 0.4, 1],
                    ["KlClusterS", 10, 0.4, 1],
                    ["KlClusterS", 15, 0.4, 1],
                    ["KlClusterS", 5, 0.5, 1],
                    ["KlClusterS", 10, 0.5, 1],
                    ["KlClusterS", 15, 0.5, 1],
                    ["KlClusterS", 5, 0.6, 1],
                    ["KlClusterS", 10, 0.6, 1],
                    ["KlClusterS", 15, 0.6, 1],
                    ["KlClusterS", 5, 0.7, 1],
                    ["KlClusterS", 10, 0.7, 1],
                    ["KlClusterS", 15, 0.7, 1],
                    ["KlClusterS", 5, 0.8, 1],
                    ["KlClusterS", 10, 0.8, 1],
                    ["KlClusterS", 15, 0.8, 1],
                    ["KlClusterS", 5, 0.9, 1],
                    ["KlClusterS", 10, 0.9, 1],
                    ["KlClusterS", 15, 0.9, 1],
                    ["KlClusterS", 5, 1.0, 1],
                    ["KlClusterS", 10, 1.0, 1],
                    ["KlClusterS", 15, 1.0, 1] ,
                    ["KlClusterS", 5, 1.1, 1],
                    ["KlClusterS", 10, 1.1, 1],
                    ["KlClusterS", 15, 1.1, 1] ,
                                        ["KlClusterS", 5, 1.2, 1],
                    ["KlClusterS", 10, 1.2, 1],
                    ["KlClusterS", 15, 1.2, 1] ,
                                        ["KlClusterS", 5, 1.3, 1],
                    ["KlClusterS", 10, 1.3, 1],
                    ["KlClusterS", 15, 1.3, 1] ,
                                        ["KlClusterS", 5, 1.4, 1],
                    ["KlClusterS", 10, 1.4, 1],
                    ["KlClusterS", 15, 1.4, 1] ,
                                        ["KlClusterS", 5, 1.5, 1],
                    ["KlClusterS", 10, 1.5, 1],
                    ["KlClusterS", 15, 1.5, 1] ,
                    ["KlClusterS", 5, 2.0, 1],
                    ["KlClusterS", 10, 2.0, 1],
                    ["KlClusterS", 15, 2.0, 1] ,
                ], 5],
                [[TAG for TAG in range(1, 15)], [
                    ["KlClusterF", 5, 1.0, 0.2],
                    ["KlClusterF", 10, 1.0, 0.2],
                    ["KlClusterF", 15, 1.0, 0.2],
                    ["KlClusterF", 5, 1.0, 0.4],
                    ["KlClusterF", 10, 1.0, 0.4],
                    ["KlClusterF", 15, 1.0, 0.4],
                    ["KlClusterF", 5, 1.0, 0.6],
                    ["KlClusterF", 10, 1.0, 0.6],
                    ["KlClusterF", 15, 1.0, 0.6],
                    ["KlClusterF", 5, 1.0, 0.8],
                    ["KlClusterF", 10, 1.0, 0.8],
                    ["KlClusterF", 15, 1.0, 0.8],
                    ["KlClusterF", 5, 1.0, 1.0],
                    ["KlClusterF", 10, 1.0, 1.0],
                    ["KlClusterF", 15, 1.0, 1.0],
                    ["KlClusterF", 5, 1.0, 1.2],
                    ["KlClusterF", 10, 1.0, 1.2],
                    ["KlClusterF", 15, 1.0, 1.2],
                    ["KlClusterF", 5, 1.0, 1.4],
                    ["KlClusterF", 10, 1.0, 1.4],
                    ["KlClusterF", 15, 1.0, 1.4],
                    ["KlClusterF", 5, 1.0, 1.6],
                    ["KlClusterF", 10, 1.0, 1.6],
                    ["KlClusterF", 15, 1.0, 1.7],
                    ["KlClusterF", 5, 1.0, 1.8],
                    ["KlClusterF", 10, 1.0, 1.8],
                    ["KlClusterF", 15, 1.0, 1.8],
                    ["KlClusterF", 5, 1.0, 2.0],
                    ["KlClusterF", 10, 1.0, 2.0],
                    ["KlClusterF", 15, 1.0, 2.0],
                    ["KlClusterF", 5, 1.0, 3.0],
                    ["KlClusterF", 10, 1.0, 3.0],
                    ["KlClusterF", 15, 1.0, 3.0],
                ], 5],

 [[TAG for TAG in range(1, 15)], [
                    ["KlClusterC", 1, 1.0, 1.0],
                    ["KlClusterC", 2, 1.0, 1.0],
                    ["KlClusterC", 3, 1.0, 1.0],
                    ["KlClusterC", 4, 1.0, 1.0],
                    ["KlClusterC", 5, 1.0, 1.0],
                    ["KlClusterC", 6, 1.0, 1.0],
                    ["KlClusterC", 7, 1.0, 1.0],
                    ["KlClusterC", 8, 1.0, 1.0],
                    ["KlClusterC", 9, 1.0, 1.0],
                    ["KlClusterC", 10, 1.0, 1.0],
                    ["KlClusterC", 11, 1.0, 1.0],
                    ["KlClusterC", 12, 1.0, 1.0],
                    ["KlClusterC", 13, 1.0, 1.0],
                    ["KlClusterC", 14, 1.0, 1.0],
                    ["KlClusterC", 15, 1.0, 1.0],
                    ["KlClusterC", 16, 1.0, 1.0],
                    ["KlClusterC", 17, 1.0, 1.0],
                    ["KlClusterC", 18, 1.0, 1.0],
                    ["KlClusterC", 19, 1.0, 1.0],
                    ["KlClusterC", 20, 1.0, 1.0],
                    ["KlClusterC", 30, 1.0, 1.0],
                ], 5],

                [[TAG for TAG in range(1, 15)], [["KlCluster", 4, 0.8, 1.25]], 5],
                [[TAG for TAG in range(1, 15)], [["KlCluster", 5, 0.8, 1.35]], 5],
                [[TAG for TAG in range(1, 15)], [["KlCluster", 6, 0.85, 1.4]], 5],
                [[TAG for TAG in range(1, 15)], [["KlCluster", 7, 0.8, 0.95]], 5],
                [[TAG for TAG in range(1, 15)], [["KlCluster", 8, 0.8, 1.4]], 3],
                [[TAG for TAG in range(1, 15)], [["KlCluster", 9, 0.8, 0.95]], 3],
                [[TAG for TAG in range(1, 15)], [["KlCluster", 10, 0.8, 1.7]], 3],
                [[TAG for TAG in range(1, 15)], [["KlCluster", 11, 0.8, 0.95]], 3],
                [[TAG for TAG in range(1, 15)], [["KlCluster", 12, 0.8, 0.8]], 3],
                [[TAG for TAG in range(1, 15)], [["KlCluster", 13, 0.8, 1.7]], 3],
                [[TAG for TAG in range(1, 15)], [["KlCluster", 14, 0.85, 0.85]], 3],
                [[TAG for TAG in range(1, 15)], [["KlCluster", 15, 0.8, 1.05]], 3],
                [[TAG for TAG in range(1, 15)], [["KlCluster", 16, 0.8, 1.1]], 3],
                [[TAG for TAG in range(1, 15)], [["KlCluster", 17, 0.8, 1.5]], 3],
                [[TAG for TAG in range(1, 15)], [["KlCluster", 18, 0.8, 1.05]], 3],
                [[TAG for TAG in range(1, 15)], [["KlCluster", 19, 0.8, 1.1]], 3],
                [[TAG for TAG in range(1, 15)], [["Aca"], ["Haca"], ["Gmm"]], 3],
                [[TAG for TAG in range(1, 15)], [["Tmm"]], 3]

                ]


for instances, methods, num_iterations in configurations:
    for TAG in instances:
            for method in methods:
                print(method)
                complexity, simpDelta, freeDelta = None, None, None
                init_time = 0
                if "KlCluster" in method[0]:
                    _, complexity, simpDelta, freeDelta = method

                    print(num_iterations)
                    print(method)
                    #KlClusterCMUSolver(TAG, simpDelta, freeDelta, complexity, num_iterations=ITERATIONS)
                    
                    solver = KlClusterCMUSolver(TAG, simpDelta, freeDelta, complexity, num_iterations=ITERATIONS)
                    segments = solver.solve()
                if method[0] == "Tmm":
                    solver = TmmCMUSolver(TAG, num_iterations=ITERATIONS)
                    segments = solver.solve()
                if method[0] == "Aca":
                    solver = AcaCMUSolver(TAG, "aca", num_iterations=ITERATIONS)
                    segments = solver.solve()
                if method[0] == "Haca":
                    solver = AcaCMUSolver(TAG, "haca", num_iterations=ITERATIONS)
                    segments = solver.solve()
                if method[0] == "Gmm":
                    solver = AcaCMUSolver(TAG, "gmm", num_iterations=ITERATIONS)
                    segments = solver.solve()

                new_entry = {
                        "TAG" : TAG, 
                        "Method" : method[0], 
                        "nClusters" : len(solver.clusters), 
                        "Time" : solver.getExecutionTime() + solver.getInitTime(),
                        "SolvingTime" : solver.getExecutionTime(),
                        "SimplifyTime" : solver.getInitTime(),
                        "Accuracy" : solver.calculateAccuracy(segments, TAG), 
                        "True Accuracy" : solver.calculateAccuracyTrueLabels(segments, TAG), 
                        "Macro Precision" : solver.calculateMacroPrecision(segments, TAG), 
                        "Macro Recall" : solver.calculateMacroRecall(segments, TAG), 
                        "Macro F1" : solver.calculateMacroF1(segments, TAG), 
                        "Num iterations" : num_iterations,
                        "Complexity" : complexity,
                        "Simp Delta" : simpDelta,
                        "Free Delta" : freeDelta,
                        #"OD" : solver.getOD() if "KlCluster" in method[0] else 0,
                    }
                new_entry_df = pd.DataFrame([new_entry])

                df = pd.concat([df, new_entry_df], ignore_index=True)
                df.to_csv(csv_file, index=False)

