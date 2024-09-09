from cmuSolvers import KlClusterCMUSolver
import numpy as np
import sys
import os
import json
from random import random
import time
from multiprocessing import Process, Lock

SIMP_DELTAS = np.arange(0.8, 2, 0.05)
FREE_DELTAS = np.arange(0.8, 2, 0.05)
COMPLEXITIES = np.arange(4, 20, 1)
# SIMP_DELTAS = [1.25]
# FREE_DELTAS = [1.25]
# COMPLEXITIES = [10]

def main():
    starttime = time.time()
    resultsfilepath = os.path.join(os.path.dirname(__file__), "gridsearch.csv")

    print(f"testing {len(SIMP_DELTAS) * len(FREE_DELTAS) * len(COMPLEXITIES)} configurations")

    nthreads = 6
    print(f"Using {nthreads} threads")

    with open(resultsfilepath, "w") as f:
        f.write("TAG,SIMP_DELTA,FREE_DELTA,COMPLEXITY,acc,accTrue,macroPrec,macroRec,macroF1,nClusters\n")

    # lock to sync file writing
    filelock = Lock()

    # start processes
    processes = [Process(target=task, args=(nthreads, i, resultsfilepath, filelock)) for i in range(nthreads)]
    for p in processes:
        p.start()

    # wait for all processes to finish
    for p in processes:
        p.join()

    print(f"Done, took {time.time() - starttime} seconds.")
 
def task(nthreads, threadidx, resultsfilepath, filelock):
    print(f"Started thread {threadidx}/{nthreads}")
    counter = 0
    for s_delta in SIMP_DELTAS:
        for f_delta in FREE_DELTAS:
            for l in COMPLEXITIES:
                counter += 1
                if s_delta>f_delta:
                    continue
                if counter % nthreads != threadidx:
                    continue
                
                # run experiment in separate process to catch assertion errors
                p = Process(target=runsingleconfiguration, args=(s_delta, f_delta, l, resultsfilepath, filelock))
                p.start()
                p.join()

                if counter % 10 == 0:
                    print(f"Thread {threadidx} Progress: {counter/(len(SIMP_DELTAS) * len(FREE_DELTAS) * len(COMPLEXITIES))*100}%")

    print(f"Thread {threadidx}/{nthreads} finished")

def runsingleconfiguration(s_delta, f_delta, l, resultsfilepath, filelock):
    result = {
        "nClusters": [],
        "acc": [],
        "accTrue": [],
        "macroPrec": [],
        "macroRec": [],
        "macroF1": []
    }

    for TAG in range(1, 15):
        print(f"TAG {TAG}")
        print(f"{s_delta} {f_delta} {l}")
        sys.stdout.flush() # flush stdout of subprocess
        sys.stderr.flush()

        solver = KlClusterCMUSolver(TAG, SIMP_DELTA=s_delta, FREE_DELTA=f_delta, COMPLEXITY=l)
        segments = solver.solve(mergeOverlappingClusters=False)

        acc = solver.calculateAccuracy(segments, TAG)
        accTrue = solver.calculateAccuracyTrueLabels(segments, TAG)
        macroPrec = solver.calculateMacroPrecision(segments, TAG)
        macroRec = solver.calculateMacroRecall(segments, TAG)
        macroF1 = solver.calculateMacroF1(segments, TAG)

        nClusters = len(solver.clusters)
        
        with filelock:
            with open(resultsfilepath, "a") as f:
                result_json = json.dumps(result, separators=(',', ':'))
                result_json.replace("\"", "\"\"")
                f.write(f"{TAG},{round(s_delta, 5)},{round(f_delta, 5)},{l},{acc},{accTrue},{macroPrec},{macroRec},{macroF1},{nClusters}\n")

if __name__ == '__main__':
    main()