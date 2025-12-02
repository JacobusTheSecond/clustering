import math

def refine_solution(elements, Delta, c, verbose=4):
    """
    Greedy refinement of a set cover solution.

    Verbose levels:
        5 = internal (never printed for user) – used for per-element testing
        4 = full debug (intervals, element tests, cost details)
        3 = iteration summaries + merged intervals + cost
        2 = iteration summaries + cost only
        1 = only final result
        0 = silent
    """

    # ------------------------------------------------------------
    # Print helper
    # ------------------------------------------------------------
    def vprint(level, *args, **kwargs):
        if verbose >= level:
            print(*args, **kwargs)

    # ------------------------------------------------------------
    # Start message
    # ------------------------------------------------------------
    vprint(4, "\n==========================================================")
    vprint(4, "  STARTING REFINEMENT PROCEDURE")
    vprint(4, "==========================================================\n")

    # ------------------------------------------------------------
    # 1. Infer N and M
    # ------------------------------------------------------------
    vprint(4, "==========================================================")
    vprint(4, " INFERRING N AND M FROM ELEMENTS")
    vprint(4, "==========================================================")

    all_ns = set()
    max_b = {}

    for eid, element in enumerate(elements):
        vprint(4, f"[DEBUG] Scanning element {eid}: {element}")
        for (n, a, b) in element:
            all_ns.add(n)
            max_b[n] = max(max_b.get(n, b), b)

    N = max(all_ns) + 1
    M = [max_b.get(i, 0.0) for i in range(N)]

    vprint(4, "\n[RESULT] Inferred:")
    vprint(4, f"  N = {N}")
    vprint(4, f"  M = {M}")
    vprint(4, "==========================================================\n")

    # ------------------------------------------------------------
    # 2. Weight function
    # ------------------------------------------------------------
    def f(i, a, b):
        lo = math.ceil(a)
        hi = math.floor(b)
        if hi < lo:
            return 0
        return hi - lo + 1

    # ------------------------------------------------------------
    # 3. Initialization
    # ------------------------------------------------------------
    S = set()
    covered_intervals = {i: [] for i in range(N)}
    W_total = sum(f(i, 0, M[i]) for i in range(N))

    vprint(4, "==========================================================")
    vprint(4, " INITIALIZATION")
    vprint(4, "==========================================================")
    vprint(4, f"Total weight W_total = {W_total}")
    for eid, e in enumerate(elements):
        vprint(4, f"Element {eid}: {e}, Δ = {Delta[eid]}")
    vprint(4, "==========================================================\n")

    # ------------------------------------------------------------
    # Helper: merge intervals
    # ------------------------------------------------------------
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

    # ------------------------------------------------------------
    # Helper: uncovered weight
    # ------------------------------------------------------------
    def compute_uncovered_weight():
        if verbose >= 4:
            vprint(4, "\n[DEBUG] Computing uncovered weights...")

        total_uncovered = 0
        for i in range(N):
            merged = merge_intervals(covered_intervals[i])
            weight_cov = sum(f(i, a, b) for (a, b) in merged)
            total = f(i, 0, M[i])
            uncovered = total - weight_cov

            if verbose >= 4:
                vprint(4, f"Interval {i}: merged={merged}, covered={weight_cov}, total={total}, uncovered={uncovered}")
            elif verbose == 3:
                vprint(3, f" Interval {i}: covered={weight_cov}/{total}")

            total_uncovered += uncovered

        return total_uncovered

    # ------------------------------------------------------------
    # Helper: cost
    # ------------------------------------------------------------
    def compute_cost():
        max_delta = max((Delta[e] for e in S), default=0)

        if verbose >= 4:
            vprint(4, f"[DEBUG] Computing cost, max Δ = {max_delta}")

        uncovered = compute_uncovered_weight()
        normalized = uncovered / W_total
        cost_val = c[0] * len(S) + c[1] * max_delta + c[2] * normalized

        if verbose >= 4:
            vprint(4, f"[DEBUG] Cost = {cost_val}")

        return cost_val

    # ------------------------------------------------------------
    # 4. Greedy refinement
    # ------------------------------------------------------------
    iteration = 0
    improved = True

    while improved:
        iteration += 1
        if verbose >= 2:
            print(f"\n---- ITERATION {iteration} ----")

        current_cost = compute_cost()
        if verbose >= 2:
            print(f"Current cost = {current_cost}")

        improved = False
        best_cost = current_cost
        best_element = None

        # Try each unused element
        for e_id, e_intervals in enumerate(elements):
            if e_id in S:
                continue

            # Tentatively add
            for (n, a, b) in e_intervals:
                covered_intervals[n].append((a, b))
            S.add(e_id)

            test_cost = compute_cost()
            vprint(5, f"[TEST] cost if add {e_id} = {test_cost}")

            if test_cost < best_cost:
                best_cost = test_cost
                best_element = e_id

            # Undo
            S.remove(e_id)
            for (n, a, b) in e_intervals:
                covered_intervals[n].pop()

        if best_element is not None:
            # Commit improvement
            S.add(best_element)
            for (n, a, b) in elements[best_element]:
                covered_intervals[n].append((a, b))

            improved = True

            if verbose >= 2:
                print(f"Improved → new cost = {best_cost}")

        else:
            if verbose >= 2:
                print(f"No improvement in iteration {iteration}. Stopping.")
            break

    # ------------------------------------------------------------
    # Final result
    # ------------------------------------------------------------
    final_cost = compute_cost()

    if verbose >= 1:
        print("\n==== FINAL RESULT ====")
        print(f"Selected set S = {sorted(S)}")
        print(f"Final cost     = {final_cost}")

    return S, final_cost

