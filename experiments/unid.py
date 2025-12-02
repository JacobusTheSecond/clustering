import os
import multiprocessing
import psutil

import refiner
import refiner2
import refiner3

klcluster = True
datafolder = "/home/jacobus/data/unid"
drifterFiles = [os.path.join(datafolder, file) for file in os.listdir(datafolder) if file.endswith(".txt")]

resultdict = {}

# --- Subprocess function ---
def run_solver_single_threaded(drifterFiles, delta, c2, l):
    # Force single-threaded OpenMP/BLAS
    os.environ["OMP_NUM_THREADS"] = "1"
    os.environ["OPENBLAS_NUM_THREADS"] = "1"
    os.environ["MKL_NUM_THREADS"] = "1"
    os.environ["NUMEXPR_NUM_THREADS"] = "1"
    os.environ["VECLIB_MAXIMUM_THREADS"] = "1"
    os.environ["BLIS_NUM_THREADS"] = "1"

    # Import solver inside subprocess to respect env vars
    from drifterSolvers import KlClusterDriftersSolver
    import refiner3
    import time

    # Construct the solver (not timed)
    solver = KlClusterDriftersSolver(
        drifterFiles,
        freedelta=delta,
        simpdelta=delta / 8,
        withaggressive=False,
        rossbyRadius=False,
        complexity=l,
    ) if klcluster else None

    process = psutil.Process(os.getpid())

    # --- Start timing memory and runtime for solution + refinement ---
    start = time.time()
    mem_before = process.memory_info().rss

    # Compute solution size
    size = solver.solutionSize(False)

    # Build elements
    elements = []
    for c in solver.clusters:
        element = []
        for m in c:
            element.append([
                m.curve,
                solver.cc.mapToBase(m.curve, m.start).value,
                solver.cc.mapToBase(m.curve, m.end).value
            ])
        elements.append(element)

    # Refine
    S, cost, _ = refiner3.refine_cut_delta(elements, [delta] * len(elements), (1, c2, 362), 2)

    mem_after = process.memory_info().rss
    end = time.time()

    mem_usage_mb = (mem_after - mem_before) / (1024 ** 2)
    elapsed_time = end - start

    return size, len(S), cost, mem_usage_mb, elapsed_time

# --- Main loop ---
if __name__ == "__main__":
    ctx = multiprocessing.get_context("spawn")  # fresh process to respect single-threading
    cDict = {}
    for c in [3,30,300]:

        l = 2
        resultdict = {}
        for j in range(6):
            delta = 32
            lDict = {}
            for i in range(6):
                # Run solver in single-threaded subprocess
                with ctx.Pool(processes=1) as pool:
                    size, filteredSize, cost, mem_usage, elapsed_time = pool.apply(run_solver_single_threaded, (drifterFiles, delta, c, l))

                lDict[delta] = {
                    "size": size,
                    "filteredSize": filteredSize,
                    "cost": cost,   
                    "mem_MB": mem_usage,
                    "time_s": elapsed_time
                }

                delta /= 2
                print(lDict)

            resultdict[l] = lDict
            l *= 2
            print(resultdict)
        cDict[c] = resultdict
    print("=============")
    print(cDict)
