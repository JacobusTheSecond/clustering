
import pandas as pd
import os
import timeit
import json


result = {}

ITERATIONS = 1

columns = ["points", "complexity", "SimpDelta", "FreeDelta", "Time", "SolvingTime", "SimpliTime", "NumClusters", "NumCurves", "allClusters", "OD"]

filename = 'drifters.csv'

if os.path.exists(filename) and os.path.getsize(filename) > 0:
    df = pd.read_csv(filename)

    if df.empty:
        df = pd.DataFrame(columns=columns)
else:
    df = pd.DataFrame(columns=columns)

from drifterSolvers import KlClusterDriftersSolver

n = 2000

klcluster = True
datafolder = os.path.join(os.path.dirname(__file__), "../data_drifters/world3d_txt")
drifterFiles = [os.path.join(datafolder, file) for file in os.listdir(datafolder) if file.endswith(".txt")]


configurations = [

                [3, 30_000, 80_000],
                [3, 30_000, 160_000],
                [3, 30_000, 240_000],
                [3, 30_000, 320_000],
                
                [3, 30_000, 80_000],
                [3, 60_000, 160_000],
                [3, 90_000, 240_000],
                [3, 120_000, 320_000],

                [3, 30_000, 160_000],
                [3, 60_000, 160_000],
                [3, 90_000, 160_000],
                [3, 120_000, 160_000],

                [1, 60_000, 160_000],
                [3, 60_000, 160_000],
                [5, 60_000, 160_000],
                [10, 60_000, 160_000],
                # [[TAG for TAG in range(1, 15)], [["KlCluster", 5, 0.8, 1.35]], 5],
                # [[TAG for TAG in range(1, 15)], [["KlCluster", 6, 0.85, 1.4]], 5],
                # [[TAG for TAG in range(1, 15)], [["KlCluster", 7, 0.8, 0.95]], 5],
                # [[TAG for TAG in range(1, 15)], [["KlCluster", 8, 0.8, 1.4]], 3],
                # [[TAG for TAG in range(1, 15)], [["KlCluster", 9, 0.8, 0.95]], 3],
                # [[TAG for TAG in range(1, 15)], [["KlCluster", 10, 0.8, 1.7]], 3],
                # [[TAG for TAG in range(1, 15)], [["KlCluster", 11, 0.8, 0.95]], 3],
                # [[TAG for TAG in range(1, 15)], [["KlCluster", 12, 0.8, 0.8]], 3],
                # [[TAG for TAG in range(1, 15)], [["KlCluster", 13, 0.8, 1.7]], 3],
                # [[TAG for TAG in range(1, 15)], [["KlCluster", 14, 0.85, 0.85]], 3],
                # [[TAG for TAG in range(1, 15)], [["KlCluster", 15, 0.8, 1.05]], 3],
                # [[TAG for TAG in range(1, 15)], [["KlCluster", 16, 0.8, 1.1]], 3],
                # [[TAG for TAG in range(1, 15)], [["KlCluster", 17, 0.8, 1.5]], 3],
                # [[TAG for TAG in range(1, 15)], [["KlCluster", 18, 0.8, 1.05]], 3],
                # [[TAG for TAG in range(1, 15)], [["KlCluster", 19, 0.8, 1.1]], 3],
                # [[TAG for TAG in range(1, 15)], [["Aca"], ["Haca"], ["Gmm"]], 3],
                # [[TAG for TAG in range(1, 15)], [["Tmm"]], 3]

                ]

times = [1e4, 2e4, 3e4, 4e4, 5e4, 6e4, 7e4, 8e4, 9e4, 
        1e5, 2e5, 3e5, 4e5, 5e5, 6e5, 7e5, 8e5, 9e5, 
        1e6, 2e6, 3e6, 4e6, 5e6, 6e6, 7e6, 8e6, 9e6, 
        1e7, 2e7]


for complexity, simpDELTA, freeDELTA in configurations:
    for allClusters in [0, 1]:    
        for points in [5e5]:#times:
            print(points)
            #solver = KlClusterCMUSolver(TAG, simpDELTA, freeDELTA, complexity)
            #init_time = (timeit.timeit(lambda: KlClusterCMUSolver(TAG, simpDelta, freeDelta, complexity), number=num_iterations) / num_iterations) if num_iterations != 0 else 0
                            
            solver = KlClusterDriftersSolver(drifterFiles, points ,simpDELTA, freeDELTA, complexity, iterations=ITERATIONS, storeFiles=True)
            solver.solve(allClusters)
            
            print("f")
            new_entry = {"points": int((points//1000)*1000), 
        "complexity" : complexity, 
        "SimpDelta" : simpDELTA, 
        "FreeDelta" : freeDELTA, 
        "Time" : solver.getExecutionTime() + solver.getInitTime(),
        "SolvingTime" : solver.getExecutionTime(), 
        "SimpliTime" : solver.getInitTime(), 
        "NumClusters" : len(solver.clustercurves), 
        "NumCurves" : len(solver.datacurves), 
        "allClusters": allClusters,
        "OD": solver.getOD()}

            
            #if solver.getExecutionTime() + solver.getInitTime():

            new_entry_df = pd.DataFrame([new_entry])
            print(new_entry_df)

            df = pd.concat([df, new_entry_df], ignore_index=True)
            df.to_csv(filename, index=False)
            if solver.getExecutionTime() + solver.getInitTime() > 300:
                print("1111")
                solver.plotInput()
                solver.plotResults(relevant=True)
                solver.plotResults(relevant=False)
                solver.plotInputAndResult(relevant=True)
                solver.plotInputAndResult(relevant=False)
                solver.plotCurve(relevant=False)
                solver.plotCurve(relevant=True)
                solver.drawInputField()
                solver.drawResultField(relevant=True)
                solver.drawResultField(relevant=False)

