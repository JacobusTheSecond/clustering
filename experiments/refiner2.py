import math
import random

###############################################################################
# Utility: Weight function
###############################################################################

def f(i, a, b):
    """Weight = number of integers in [a, b]."""
    lo = math.ceil(a)
    hi = math.floor(b)
    if hi < lo:
        return 0
    return hi - lo + 1


###############################################################################
# Merge intervals (old approach)
###############################################################################

def merge_intervals(intervals):
    if not intervals:
        return []
    intervals = sorted(intervals, key=lambda x: x[0])
    merged = [intervals[0]]
    for a, b in intervals[1:]:
        la, lb = merged[-1]
        if a <= lb:
            merged[-1] = (la, max(lb, b))
        else:
            merged.append((a, b))
    return merged


###############################################################################
# Cut intervals to maintain disjoint coverage (new approach)
###############################################################################

def cut_interval(existing, a, b):
    """Returns list of subintervals of (a,b) not covered by existing."""
    result = [(a, b)]
    for (c, d) in existing:
        new_result = []
        for (x, y) in result:
            if y <= c or x >= d:
                new_result.append((x, y))
                continue
            if x < c:
                new_result.append((x, c))
            if y > d:
                new_result.append((d, y))
        result = new_result
        if not result:
            break
    return result


###############################################################################
# Print helper
###############################################################################

def vprint(level, verbose, *args, **kwargs):
    if verbose >= level:
        print(*args, **kwargs)


###############################################################################
# Infer N and M automatically from elements
###############################################################################

def infer_N_M(elements):
    all_ns = set()
    max_b = {}
    for element in elements:
        for (n, a, b) in element:
            all_ns.add(n)
            max_b[n] = max(max_b.get(n, b), b)
    N = max(all_ns) + 1
    M = [max_b.get(i, 0.0) for i in range(N)]
    return N, M


###############################################################################
# Merge-based refinement
###############################################################################

def refine_merge(elements, Delta, c, verbose=4):
    N, M = infer_N_M(elements)
    S = set()
    covered = {i: [] for i in range(N)}

    iteration = 0
    improved = True

    while improved:
        iteration += 1
        vprint(3, verbose, f"\n--- Merge-based iteration {iteration} ---")
        current_cost = compute_cost_merge(S, covered, M, Delta, c)
        vprint(2, verbose, f"Current cost = {current_cost}")

        improved = False
        best_cost = current_cost
        best_e = None

        for e_id, e_list in enumerate(elements):
            if e_id in S:
                continue

            # Tentatively add
            for (n, a, b) in e_list:
                covered[n].append((a, b))

            test_cost = compute_cost_merge(S | {e_id}, covered, M, Delta, c)

            if verbose >= 4:
                vprint(4, verbose, f"Test element {e_id}, cost = {test_cost}, intervals = {e_list}")

            # Undo
            for (n, a, b) in e_list:
                covered[n].pop()

            if test_cost < best_cost:
                best_cost = test_cost
                best_e = e_id

        if best_e is not None:
            S.add(best_e)
            for (n, a, b) in elements[best_e]:
                covered[n].append((a, b))
            improved = True
            vprint(2, verbose, f"Iteration {iteration} improved → new cost = {best_cost}")
            if verbose >= 3:
                for i in range(N):
                    merged = merge_intervals(covered[i])
                    vprint(3, verbose, f"Interval {i} coverage: {merged}")
        else:
            vprint(2, verbose, f"No improvement in iteration {iteration}. Stopping.")
            break

    final_cost = compute_cost_merge(S, covered, M, Delta, c)
    vprint(1, verbose, f"\nMerge-based final S = {sorted(S)}, cost = {final_cost}")
    return S, final_cost, covered


def compute_cost_merge(S, covered, M, Delta, c):
    max_delta = max((Delta[e] for e in S), default=0)
    W_total = sum(f(i, 0, M[i]) for i in range(len(M)))
    covered_weight = 0
    for i in range(len(M)):
        merged = merge_intervals(covered[i])
        covered_weight += sum(f(i, a, b) for (a, b) in merged)
    uncovered = W_total - covered_weight
    normalized = uncovered / W_total
    return c[0]*len(S) + c[1]*max_delta + c[2]*normalized


###############################################################################
# Cut-based refinement
###############################################################################

def refine_cut(elements, Delta, c, verbose=4):
    N, M = infer_N_M(elements)
    S = set()
    covered = {i: [] for i in range(N)}

    iteration = 0
    improved = True

    while improved:
        iteration += 1
        vprint(3, verbose, f"\n--- Cut-based iteration {iteration} ---")
        current_cost = compute_cost_cut(S, covered, M, Delta, c)
        vprint(2, verbose, f"Current cost = {current_cost}")

        improved = False
        best_cost = current_cost
        best_e = None
        best_cut_pieces = None

        for e_id, e_list in enumerate(elements):
            if e_id in S:
                continue

            cut_pieces = []
            for (n, a, b) in e_list:
                new_parts = cut_interval(covered[n], a, b)
                for (u, v) in new_parts:
                    cut_pieces.append((n, u, v))

            test_covered_weight = sum(f(i, a, b) for i in range(N) for (a, b) in covered[i]) \
                                  + sum(f(n, a, b) for (n, a, b) in cut_pieces)
            W_total = sum(f(i, 0, M[i]) for i in range(N))
            uncovered = W_total - test_covered_weight
            normalized = uncovered / W_total
            test_cost = c[0]*(len(S)+1) + c[1]*max(Delta[e] for e in (S|{e_id})) + c[2]*normalized

            if verbose >= 4:
                vprint(4, verbose, f"Test element {e_id}, cut pieces = {cut_pieces}, test_cost = {test_cost}")

            if test_cost < best_cost:
                best_cost = test_cost
                best_e = e_id
                best_cut_pieces = cut_pieces

        if best_e is not None:
            S.add(best_e)
            for (n, a, b) in best_cut_pieces:
                covered[n].append((a, b))
            improved = True
            vprint(2, verbose, f"Iteration {iteration} improved → new cost = {best_cost}")
            if verbose >= 3:
                for i in range(N):
                    vprint(3, verbose, f"Interval {i} coverage: {covered[i]}")
        else:
            vprint(2, verbose, f"No improvement in iteration {iteration}. Stopping.")
            break

    final_cost = compute_cost_cut(S, covered, M, Delta, c)
    vprint(1, verbose, f"\nCut-based final S = {sorted(S)}, cost = {final_cost}")
    return S, final_cost, covered


def compute_cost_cut(S, covered, M, Delta, c):
    max_delta = max((Delta[e] for e in S), default=0)
    W_total = sum(f(i, 0, M[i]) for i in range(len(M)))
    covered_weight = sum(f(i, a, b) for i in range(len(M)) for (a, b) in covered[i])
    uncovered = W_total - covered_weight
    normalized = uncovered / W_total
    return c[0]*len(S) + c[1]*max_delta + c[2]*normalized

