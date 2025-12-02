# ============================================================
# Standalone LaTeX Table Generator (No Input Files Needed)
# ============================================================

# ----------------------------
# ✅ INPUT FORMAT:
# DATA[c2][l][delta] = {
#   size, filteredSize, cost, mem_MB, time_s
# }
# ----------------------------

DATA = {10: {2: {16: {'size': 107, 'filteredSize': 73, 'cost': 241.3780398319049, 'mem_MB': 60904.953125, 'time_s': 207.6700439453125}, 8.0: {'size': 317, 'filteredSize': 176, 'cost': 295.2867900502265, 'mem_MB': 32237.07421875, 'time_s': 155.82630109786987}, 4.0: {'size': 893, 'filteredSize': 0, 'cost': 2011.0, 'mem_MB': 23908.671875, 'time_s': 128.4687533378601}, 2.0: {'size': 2555, 'filteredSize': 0, 'cost': 2011.0, 'mem_MB': 18940.21484375, 'time_s': 182.0060579776764}, 1.0: {'size': 7479, 'filteredSize': 0, 'cost': 2011.0, 'mem_MB': 14744.89453125, 'time_s': 430.8028619289398}, 0.5: {'size': 20061, 'filteredSize': 0, 'cost': 2011.0, 'mem_MB': 11556.41015625, 'time_s': 1572.7099568843842}}, 4: {16: {'size': 106, 'filteredSize': 73, 'cost': 241.86842335809607, 'mem_MB': 60986.48828125, 'time_s': 209.11871004104614}, 8.0: {'size': 316, 'filteredSize': 176, 'cost': 294.2645031111874, 'mem_MB': 32420.484375, 'time_s': 157.2048363685608}, 4.0: {'size': 890, 'filteredSize': 0, 'cost': 2011.0, 'mem_MB': 24250.12109375, 'time_s': 131.44017624855042}, 2.0: {'size': 2496, 'filteredSize': 0, 'cost': 2011.0, 'mem_MB': 19300.0234375, 'time_s': 188.3933777809143}, 1.0: {'size': 6932, 'filteredSize': 0, 'cost': 2011.0, 'mem_MB': 15029.2265625, 'time_s': 421.7978091239929}, 0.5: {'size': 17258, 'filteredSize': 0, 'cost': 2011.0, 'mem_MB': 11768.4375, 'time_s': 1423.8544335365295}}, 8: {16: {'size': 108, 'filteredSize': 72, 'cost': 240.76742903920083, 'mem_MB': 60993.3125, 'time_s': 210.50212979316711}, 8.0: {'size': 315, 'filteredSize': 176, 'cost': 294.63257129560543, 'mem_MB': 32459.85546875, 'time_s': 157.06381750106812}, 4.0: {'size': 877, 'filteredSize': 0, 'cost': 2011.0, 'mem_MB': 24376.3359375, 'time_s': 129.79674220085144}, 2.0: {'size': 2430, 'filteredSize': 0, 'cost': 2011.0, 'mem_MB': 19447.4609375, 'time_s': 183.12227368354797}, 1.0: {'size': 6243, 'filteredSize': 0, 'cost': 2011.0, 'mem_MB': 15160.09375, 'time_s': 376.83587288856506}, 0.5: {'size': 13577, 'filteredSize': 0, 'cost': 2011.0, 'mem_MB': 11876.3828125, 'time_s': 1165.9906468391418}}, 16: {16: {'size': 109, 'filteredSize': 71, 'cost': 241.18022734463443, 'mem_MB': 60996.9765625, 'time_s': 212.78675937652588}, 8.0: {'size': 315, 'filteredSize': 176, 'cost': 293.742699130229, 'mem_MB': 32461.8515625, 'time_s': 162.47905564308167}, 4.0: {'size': 874, 'filteredSize': 0, 'cost': 2011.0, 'mem_MB': 24390.26171875, 'time_s': 137.171284198761}, 2.0: {'size': 2342, 'filteredSize': 0, 'cost': 2011.0, 'mem_MB': 19482.09765625, 'time_s': 187.07353234291077}, 1.0: {'size': 5204, 'filteredSize': 0, 'cost': 2011.0, 'mem_MB': 15209.875, 'time_s': 333.13044714927673}, 0.5: {'size': 9455, 'filteredSize': 0, 'cost': 2011.0, 'mem_MB': 11935.5546875, 'time_s': 751.7134966850281}}, 32: {16: {'size': 108, 'filteredSize': 71, 'cost': 241.08147734393705, 'mem_MB': 60996.86328125, 'time_s': 207.01689195632935}, 8.0: {'size': 315, 'filteredSize': 176, 'cost': 293.9244889042404, 'mem_MB': 32461.75390625, 'time_s': 157.2282485961914}, 4.0: {'size': 872, 'filteredSize': 0, 'cost': 2011.0, 'mem_MB': 24390.51953125, 'time_s': 133.2907497882843}, 2.0: {'size': 2285, 'filteredSize': 0, 'cost': 2011.0, 'mem_MB': 19486.5078125, 'time_s': 180.03224396705627}, 1.0: {'size': 4251, 'filteredSize': 0, 'cost': 2011.0, 'mem_MB': 15227.84765625, 'time_s': 280.18227314949036}, 0.5: {'size': 6358, 'filteredSize': 0, 'cost': 2011.0, 'mem_MB': 11976.49609375, 'time_s': 570.712480545044}}, 64: {16: {'size': 108, 'filteredSize': 71, 'cost': 241.08147734393705, 'mem_MB': 60996.8359375, 'time_s': 207.998952627182}, 8.0: {'size': 315, 'filteredSize': 176, 'cost': 293.9244889042404, 'mem_MB': 32461.8984375, 'time_s': 158.72473645210266}, 4.0: {'size': 872, 'filteredSize': 0, 'cost': 2011.0, 'mem_MB': 24390.5234375, 'time_s': 130.11958742141724}, 2.0: {'size': 2274, 'filteredSize': 0, 'cost': 2011.0, 'mem_MB': 19486.54296875, 'time_s': 180.88470602035522}, 1.0: {'size': 3950, 'filteredSize': 0, 'cost': 2011.0, 'mem_MB': 15231.66796875, 'time_s': 270.6943106651306}, 0.5: {'size': 4909, 'filteredSize': 0, 'cost': 2011.0, 'mem_MB': 12000.42578125, 'time_s': 483.0867249965668}}}} 

def make_all_rows_sorted(c2_data, c2_vec):
    rows = []

    # ✅ Sort by l → then by delta
    for l in sorted(c2_data.keys()):
        for delta in sorted(c2_data[l].keys()):
            info = c2_data[l][delta]

            name = f"basic-{l}"
            seconds = info["time_s"]
            gbytes = info["mem_MB"] / 1024
            clusters = info["filteredSize"]
            max_f = delta
            avg_f = delta
            kcenters = info["cost"]
            max_size = ""#info["size"]
            avg_size = ""#info["size"]

            row = (
                f"{name} & {seconds:.2f} & {gbytes:.2f} & {c2_vec} & "
                f"{clusters} & {max_f} & {avg_f} & {kcenters:.2f} & "
                f"{max_size} & {avg_size}\\\\"
            )
            rows.append(row)

    # ✅ Always append 2 placeholder rows
    #placeholder = (
    #    f"basic-na & -inf & nan & {c2_vec} & "
    #    f"0 & -inf & nan & 128.00 & -inf & nan\\\\"
    #)

    #rows.append(placeholder)
    #rows.append(placeholder)

    return rows

# ============================================================
# ✅ Generate FULL LaTeX table
# ============================================================

def generate_latex_table(c2_value, c2_data):
    c2_vec = f"(1.0, {c2_value}, 2011.0)"

    rows = make_all_rows_sorted(c2_data, c2_vec)

    header = r"""
\begin{tabular}{@{}lrrrrrrrrr@{}}
\toprule
Name & Seconds & GBytes & $c2$ & Clusters & Max Frechet & Avg Frechet & kCenters & Max Size & Avg Size\\
\midrule
""".strip()

    footer = r"""
\bottomrule
\end{tabular}
""".strip()

    return "\n".join([header] + rows + [footer])

# ============================================================
# ✅ MAIN: Generate One Table
# ============================================================

if __name__ == "__main__":

    # ✅ Choose which c2 table to print
    C2_VALUE_TO_USE = list(DATA.keys())[0]

    latex_table = generate_latex_table(
        C2_VALUE_TO_USE,
        DATA[C2_VALUE_TO_USE]
    )

    print("\n✅ GENERATED LATEX TABLE (ALL ENTRIES, SORTED BY l → delta):\n")
    print(latex_table)
