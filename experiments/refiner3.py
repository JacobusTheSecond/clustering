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
# Merge intervals helper
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
# Cut intervals to maintain disjoint coverage
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
# Merge-based refinement (for correctness comparison)
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
            for (n, a, b) in e_list:
                covered[n].append((a, b))
            test_cost = compute_cost_merge(S | {e_id}, covered, M, Delta, c)
            for (n, a, b) in e_list:
                covered[n].pop()
            if verbose >= 4:
                vprint(4, verbose, f"Test element {e_id}, cost = {test_cost}, intervals = {e_list}")
            if test_cost < best_cost:
                best_cost = test_cost
                best_e = e_id

        if best_e is not None:
            S.add(best_e)
            for (n, a, b) in elements[best_e]:
                covered[n].append((a, b))
                covered[n] = merge_intervals(covered[n])
            improved = True
            vprint(2, verbose, f"Iteration {iteration} improved → new cost = {best_cost}")
            if verbose >= 3:
                for i in range(N):
                    vprint(3, verbose, f"Interval {i} coverage: {covered[i]}")
        else:
            vprint(2, verbose, f"No improvement in iteration {iteration}. Stopping.")
            break

    final_cost = compute_cost_merge(S, covered, M, Delta, c)
    vprint(1, verbose, f"\nMerge-based final S = {sorted(S)}, cost = {final_cost}")
    return S, final_cost, covered

def compute_cost_merge(S, covered, M, Delta, c):
    max_delta = max((Delta[e] for e in S), default=0)
    W_total = sum(f(i, 0, M[i]) for i in range(len(M)))
    covered_weight = sum(f(i, a, b) for i in range(len(M)) for (a, b) in merge_intervals(covered[i]))
    uncovered = W_total - covered_weight
    normalized = uncovered / W_total
    return c[0]*len(S) + c[1]*max_delta + c[2]*normalized

###############################################################################
# Cut-based delta-cost refinement (corrected)
###############################################################################

def refine_cut_delta(elements, Delta, c, verbose=4):
    N, M = infer_N_M(elements)
    S = set()
    covered = {i: [] for i in range(N)}
    current_max_delta = 0
    W_total = sum(f(i, 0, M[i]) for i in range(N))

    iteration = 0
    improved = True

    while improved:
        iteration += 1
        vprint(3, verbose, f"\n--- Cut-based iteration {iteration} ---")

        current_cost = c[0]*len(S) + c[1]*current_max_delta + c[2]*(1 - sum(f(i,a,b) for i in range(N) for (a,b) in covered[i])/W_total)
        vprint(2, verbose, f"Current cost = {current_cost}")

        improved = False
        best_total_cost = float('inf')
        best_e = None
        best_cut_pieces = None
        best_new_delta = None

        for e_id, e_list in enumerate(elements):
            if e_id in S:
                continue
            cut_pieces = []
            for (n, a, b) in e_list:
                new_parts = cut_interval(covered[n], a, b)
                cut_pieces.extend((n, u, v) for (u,v) in new_parts)
            # tentative new coverage for this element
            tentative_covered = {i: covered[i][:] for i in range(N)}
            for (n, a, b) in cut_pieces:
                tentative_covered[n].append((a, b))
                tentative_covered[n] = merge_intervals(tentative_covered[n])
            new_max_delta = max(current_max_delta, Delta[e_id])
            new_total_cost = c[0]*(len(S)+1) + c[1]*new_max_delta + c[2]*(1 - sum(f(i,a,b) for i in range(N) for (a,b) in tentative_covered[i])/W_total)
            if verbose >= 4:
                vprint(4, verbose, f"Test element {e_id}, cut pieces = {cut_pieces}, new_total_cost = {new_total_cost}")
            if new_total_cost < best_total_cost:
                best_total_cost = new_total_cost
                best_e = e_id
                best_cut_pieces = cut_pieces
                best_new_delta = new_max_delta

        if best_e is not None and best_total_cost < current_cost:
            S.add(best_e)
            current_max_delta = best_new_delta
            for (n, a, b) in best_cut_pieces:
                covered[n].append((a, b))
                covered[n] = merge_intervals(covered[n])
            improved = True
            vprint(2, verbose, f"Iteration {iteration} improved → new cost = {best_total_cost}")
            if verbose >= 3:
                for i in range(N):
                    vprint(3, verbose, f"Interval {i} coverage: {covered[i]}")
        else:
            improved = False
            vprint(2, verbose, f"No improvement in iteration {iteration}. Stopping.")

    final_cost = c[0]*len(S) + c[1]*current_max_delta + c[2]*(1 - sum(f(i,a,b) for i in range(N) for (a,b) in covered[i])/W_total)
    vprint(1, verbose, f"\nCut-based final S = {sorted(S)}, cost = {final_cost}")
    return S, final_cost, covered

###############################################################################
# Random instance generator
###############################################################################

def random_instance(num_elements=8, pieces_per_element=3):
    num_intervals = random.randint(3, 6)
    M = [random.randint(5, 20) for _ in range(num_intervals)]
    elements = []
    for _ in range(num_elements):
        e = []
        for _ in range(pieces_per_element):
            n = random.randint(0, num_intervals-1)
            a = random.uniform(0, M[n]-1)
            b = random.uniform(a, M[n])
            e.append((n, a, b))
        elements.append(e)
    Delta = [random.uniform(0, 10) for _ in range(num_elements)]
    return elements, Delta

###############################################################################
# Test harness
###############################################################################

def test_verbose(verbose=4):
    random.seed(0)
    elements, Delta = random_instance()

    print("\n=== Merge-based refinement ===")
    S_merge, cost_merge, _ = refine_merge(elements, Delta, (1,1,1), verbose=verbose)

    print("\n=== Cut-based delta-cost refinement ===")
    S_cut, cost_cut, _ = refine_cut_delta(elements, Delta, (1,1,1), verbose=verbose)

    print("\n===============================================")
    print(" Merge-based result:", sorted(S_merge), "cost =", cost_merge)
    print(" Cut-based result  :", sorted(S_cut), "cost =", cost_cut)
    print("===============================================")

    assert S_merge == S_cut, "Selected sets differ!"
    assert abs(cost_merge - cost_cut) < 1e-6, "Costs differ!"
    print("SUCCESS: Both approaches match exactly!")

if __name__ == "__main__":
    test_verbose(verbose=4)

