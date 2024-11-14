from cmuSolvers import KlClusterCMUSolver
from cmuSolvers import AcaCMUSolver
from cmuSolvers import TmmCMUSolver
import pandas as pd
import os
import timeit
import json

file_path = "result.json"

result = {}

ITERATIONS = 10

columns = ["TAG", "Method", "nClusters", "Time", "Accuracy", "True Accuracy", "Macro Precision", "Macro Recall", "Macro F1", "Num iterations", "Complexity", "Simp Delta", "Free Delta"]

csv_file = 'output.csv'

if os.path.exists(csv_file) and os.path.getsize(csv_file) > 0:
    df = pd.read_csv(csv_file)

    if df.empty:
        df = pd.DataFrame(columns=columns)
else:
    df = pd.DataFrame(columns=columns)


configurations = [

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
                complexity, simpDelta, freeDelta = None, None, None
                init_time = 0
                if method[0] == "KlCluster":
                    _, complexity, simpDelta, freeDelta = method
                    print(num_iterations)
                    init_time = (timeit.timeit(lambda: KlClusterCMUSolver(TAG, simpDelta, freeDelta, complexity), number=num_iterations) / num_iterations) if num_iterations != 0 else 0
                    
                    solver = KlClusterCMUSolver(TAG, simpDelta, freeDelta, complexity)
                    segments = solver.solveAndTime(timing_rounds=num_iterations)
                if method[0] == "Tmm":
                    solver = TmmCMUSolver(TAG)
                    segments = solver.solveAndTime(num_iterations)
                if method[0] == "Aca":
                    solver = AcaCMUSolver(TAG, "aca")
                    segments = solver.solveAndTime(num_iterations)
                if method[0] == "Haca":
                    solver = AcaCMUSolver(TAG, "haca")
                    segments = solver.solveAndTime(num_iterations)
                if method[0] == "Gmm":
                    solver = AcaCMUSolver(TAG, "gmm")
                    segments = solver.solveAndTime(num_iterations)

                new_entry = {
                     "TAG" : TAG, 
                     "Method" : method[0], 
                     "nClusters" : len(solver.clusters), 
                     "Time" : solver.getExecutionTime() + init_time if method[0] == "KlCluster" else 0,
                     "SolvingTime" : solver.getExecutionTime(),
                     "SimplifyTime" : init_time,
                     "Accuracy" : solver.calculateAccuracy(segments, TAG), 
                     "True Accuracy" : solver.calculateAccuracyTrueLabels(segments, TAG), 
                     "Macro Precision" : solver.calculateMacroPrecision(segments, TAG), 
                     "Macro Recall" : solver.calculateMacroRecall(segments, TAG), 
                     "Macro F1" : solver.calculateMacroF1(segments, TAG), 
                     "Num iterations" : num_iterations,
                     "Complexity" : complexity,
                     "Simp Delta" : simpDelta,
                     "Free Delta" : freeDelta
                    }
                new_entry_df = pd.DataFrame([new_entry])

                df = pd.concat([df, new_entry_df], ignore_index=True)
                df.to_csv('output.csv', index=False)

