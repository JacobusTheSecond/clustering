import re
from pathlib import Path

# ----------------------------
# Paths (relative to this script)
# ----------------------------
script_dir = Path(__file__).parent
latex_path = (script_dir /
    "../discrete-subtrajectory-clustering-test-suite/data/plots/socg_table_athens_centers.tex"
).resolve()

output_path = latex_path.with_name("socg_table_athens_centers_UPDATED.tex")

# ----------------------------
# ✅ Data format:
# data[c2][l][delta] = metrics
# ----------------------------
data = {0.003: {2: {4000: {'size': 1, 'filteredSize': 1, 'cost': 13.0, 'mem_MB': 19.65234375, 'time_s': 0.07989740371704102}, 2000.0: {'size': 2, 'filteredSize': 2, 'cost': 8.0, 'mem_MB': 233.2109375, 'time_s': 0.527998685836792}, 1000.0: {'size': 7, 'filteredSize': 4, 'cost': 7.766197183098598, 'mem_MB': 421.38671875, 'time_s': 0.968900203704834}, 500.0: {'size': 22, 'filteredSize': 13, 'cost': 16.843661971830983, 'mem_MB': 894.01953125, 'time_s': 2.6379990577697754}, 250.0: {'size': 63, 'filteredSize': 21, 'cost': 32.38661971830986, 'mem_MB': 1447.375, 'time_s': 4.768705129623413}, 125.0: {'size': 99, 'filteredSize': 25, 'cost': 42.97105216778287, 'mem_MB': 1446.49609375, 'time_s': 4.8310816287994385}}, 4: {4000: {'size': 1, 'filteredSize': 1, 'cost': 13.0, 'mem_MB': 19.65234375, 'time_s': 0.08008289337158203}, 2000.0: {'size': 2, 'filteredSize': 2, 'cost': 8.0, 'mem_MB': 233.2265625, 'time_s': 0.5339741706848145}, 1000.0: {'size': 7, 'filteredSize': 4, 'cost': 7.766197183098598, 'mem_MB': 422.0078125, 'time_s': 0.953587532043457}, 500.0: {'size': 23, 'filteredSize': 12, 'cost': 15.978873239436624, 'mem_MB': 908.14453125, 'time_s': 2.672074556350708}, 250.0: {'size': 55, 'filteredSize': 21, 'cost': 30.673943661971833, 'mem_MB': 1468.5, 'time_s': 4.582108974456787}, 125.0: {'size': 89, 'filteredSize': 24, 'cost': 42.15152449770885, 'mem_MB': 1469.2890625, 'time_s': 4.897400140762329}}, 8: {4000: {'size': 1, 'filteredSize': 1, 'cost': 13.0, 'mem_MB': 19.65234375, 'time_s': 0.07731819152832031}, 2000.0: {'size': 2, 'filteredSize': 2, 'cost': 8.0, 'mem_MB': 233.203125, 'time_s': 0.533336877822876}, 1000.0: {'size': 7, 'filteredSize': 4, 'cost': 7.766197183098598, 'mem_MB': 421.765625, 'time_s': 0.9601240158081055}, 500.0: {'size': 25, 'filteredSize': 11, 'cost': 15.790140845070425, 'mem_MB': 910.43359375, 'time_s': 2.549506664276123}, 250.0: {'size': 54, 'filteredSize': 21, 'cost': 30.58380281690141, 'mem_MB': 1480.671875, 'time_s': 4.752018451690674}, 125.0: {'size': 87, 'filteredSize': 22, 'cost': 39.42963517800493, 'mem_MB': 1482.76171875, 'time_s': 4.982973098754883}}, 16: {4000: {'size': 1, 'filteredSize': 1, 'cost': 13.0, 'mem_MB': 19.65234375, 'time_s': 0.07825589179992676}, 2000.0: {'size': 2, 'filteredSize': 2, 'cost': 8.0, 'mem_MB': 233.19921875, 'time_s': 0.5174915790557861}, 1000.0: {'size': 7, 'filteredSize': 4, 'cost': 7.766197183098598, 'mem_MB': 422.203125, 'time_s': 0.9779973030090332}, 500.0: {'size': 25, 'filteredSize': 11, 'cost': 15.790140845070425, 'mem_MB': 910.48828125, 'time_s': 2.6438801288604736}, 250.0: {'size': 57, 'filteredSize': 20, 'cost': 29.854225352112678, 'mem_MB': 1480.67578125, 'time_s': 4.973145961761475}, 125.0: {'size': 84, 'filteredSize': 23, 'cost': 38.534675713782164, 'mem_MB': 1484.61328125, 'time_s': 4.94984769821167}}, 32: {4000: {'size': 1, 'filteredSize': 1, 'cost': 13.0, 'mem_MB': 19.65625, 'time_s': 0.0793614387512207}, 2000.0: {'size': 2, 'filteredSize': 2, 'cost': 8.0, 'mem_MB': 233.23828125, 'time_s': 0.537508487701416}, 1000.0: {'size': 7, 'filteredSize': 4, 'cost': 7.766197183098598, 'mem_MB': 422.74609375, 'time_s': 0.9608972072601318}, 500.0: {'size': 25, 'filteredSize': 11, 'cost': 15.790140845070425, 'mem_MB': 910.48828125, 'time_s': 2.5970544815063477}, 250.0: {'size': 57, 'filteredSize': 18, 'cost': 30.10774647887324, 'mem_MB': 1482.32421875, 'time_s': 4.87201452255249}, 125.0: {'size': 84, 'filteredSize': 23, 'cost': 38.534675713782164, 'mem_MB': 1485.33203125, 'time_s': 4.947934865951538}}, 64: {4000: {'size': 1, 'filteredSize': 1, 'cost': 13.0, 'mem_MB': 19.6484375, 'time_s': 0.07876849174499512}, 2000.0: {'size': 2, 'filteredSize': 2, 'cost': 8.0, 'mem_MB': 233.2265625, 'time_s': 0.5359518527984619}, 1000.0: {'size': 7, 'filteredSize': 4, 'cost': 7.766197183098598, 'mem_MB': 421.88671875, 'time_s': 0.980729341506958}, 500.0: {'size': 25, 'filteredSize': 11, 'cost': 15.790140845070425, 'mem_MB': 911.19921875, 'time_s': 2.5807430744171143}, 250.0: {'size': 57, 'filteredSize': 18, 'cost': 30.10774647887324, 'mem_MB': 1482.32421875, 'time_s': 4.8601768016815186}, 125.0: {'size': 84, 'filteredSize': 23, 'cost': 38.534675713782164, 'mem_MB': 1484.1640625, 'time_s': 4.937316179275513}}}, 0.03: {2: {4000: {'size': 1, 'filteredSize': 1, 'cost': 121.0, 'mem_MB': 19.625, 'time_s': 0.07735848426818848}, 2000.0: {'size': 2, 'filteredSize': 2, 'cost': 62.0, 'mem_MB': 233.19921875, 'time_s': 0.5336275100708008}, 1000.0: {'size': 7, 'filteredSize': 4, 'cost': 34.7661971830986, 'mem_MB': 420.92578125, 'time_s': 0.9616093635559082}, 500.0: {'size': 22, 'filteredSize': 13, 'cost': 30.343661971830983, 'mem_MB': 894.02734375, 'time_s': 2.670922040939331}, 250.0: {'size': 63, 'filteredSize': 21, 'cost': 39.13661971830986, 'mem_MB': 1447.015625, 'time_s': 4.776504993438721}, 125.0: {'size': 99, 'filteredSize': 25, 'cost': 46.34605216778287, 'mem_MB': 1448.0859375, 'time_s': 4.769676685333252}}, 4: {4000: {'size': 1, 'filteredSize': 1, 'cost': 121.0, 'mem_MB': 19.65234375, 'time_s': 0.07737517356872559}, 2000.0: {'size': 2, 'filteredSize': 2, 'cost': 62.0, 'mem_MB': 233.203125, 'time_s': 0.534482479095459}, 1000.0: {'size': 7, 'filteredSize': 4, 'cost': 34.7661971830986, 'mem_MB': 421.5703125, 'time_s': 0.9617054462432861}, 500.0: {'size': 23, 'filteredSize': 12, 'cost': 29.478873239436624, 'mem_MB': 908.7890625, 'time_s': 2.646043062210083}, 250.0: {'size': 55, 'filteredSize': 21, 'cost': 37.42394366197183, 'mem_MB': 1469.265625, 'time_s': 4.706991910934448}, 125.0: {'size': 89, 'filteredSize': 24, 'cost': 45.52652449770885, 'mem_MB': 1470.16796875, 'time_s': 4.895172357559204}}, 8: {4000: {'size': 1, 'filteredSize': 1, 'cost': 121.0, 'mem_MB': 19.63671875, 'time_s': 0.0784766674041748}, 2000.0: {'size': 2, 'filteredSize': 2, 'cost': 62.0, 'mem_MB': 233.19921875, 'time_s': 0.5264379978179932}, 1000.0: {'size': 7, 'filteredSize': 4, 'cost': 34.7661971830986, 'mem_MB': 422.69140625, 'time_s': 0.9619870185852051}, 500.0: {'size': 25, 'filteredSize': 11, 'cost': 29.290140845070425, 'mem_MB': 910.9296875, 'time_s': 2.6080470085144043}, 250.0: {'size': 54, 'filteredSize': 21, 'cost': 37.33380281690141, 'mem_MB': 1480.2421875, 'time_s': 4.75508713722229}, 125.0: {'size': 87, 'filteredSize': 22, 'cost': 42.80463517800493, 'mem_MB': 1482.0546875, 'time_s': 4.9571692943573}}, 16: {4000: {'size': 1, 'filteredSize': 1, 'cost': 121.0, 'mem_MB': 19.640625, 'time_s': 0.07775712013244629}, 2000.0: {'size': 2, 'filteredSize': 2, 'cost': 62.0, 'mem_MB': 233.2109375, 'time_s': 0.5278255939483643}, 1000.0: {'size': 7, 'filteredSize': 4, 'cost': 34.7661971830986, 'mem_MB': 421.51953125, 'time_s': 0.9655723571777344}, 500.0: {'size': 25, 'filteredSize': 11, 'cost': 29.290140845070425, 'mem_MB': 910.51953125, 'time_s': 2.593430519104004}, 250.0: {'size': 57, 'filteredSize': 20, 'cost': 36.60422535211268, 'mem_MB': 1480.71875, 'time_s': 4.911503076553345}, 125.0: {'size': 84, 'filteredSize': 23, 'cost': 41.909675713782164, 'mem_MB': 1485.6328125, 'time_s': 4.969274520874023}}, 32: {4000: {'size': 1, 'filteredSize': 1, 'cost': 121.0, 'mem_MB': 19.62890625, 'time_s': 0.0803227424621582}, 2000.0: {'size': 2, 'filteredSize': 2, 'cost': 62.0, 'mem_MB': 233.2109375, 'time_s': 0.5319926738739014}, 1000.0: {'size': 7, 'filteredSize': 4, 'cost': 34.7661971830986, 'mem_MB': 421.5078125, 'time_s': 0.965491771697998}, 500.0: {'size': 25, 'filteredSize': 11, 'cost': 29.290140845070425, 'mem_MB': 911.33203125, 'time_s': 2.6546170711517334}, 250.0: {'size': 57, 'filteredSize': 18, 'cost': 36.85774647887324, 'mem_MB': 1482.27734375, 'time_s': 4.855437994003296}, 125.0: {'size': 84, 'filteredSize': 23, 'cost': 41.909675713782164, 'mem_MB': 1484.5390625, 'time_s': 5.010126113891602}}, 64: {4000: {'size': 1, 'filteredSize': 1, 'cost': 121.0, 'mem_MB': 19.65234375, 'time_s': 0.07872605323791504}, 2000.0: {'size': 2, 'filteredSize': 2, 'cost': 62.0, 'mem_MB': 233.19921875, 'time_s': 0.5409114360809326}, 1000.0: {'size': 7, 'filteredSize': 4, 'cost': 34.7661971830986, 'mem_MB': 421.90234375, 'time_s': 0.9660120010375977}, 500.0: {'size': 25, 'filteredSize': 11, 'cost': 29.290140845070425, 'mem_MB': 910.515625, 'time_s': 2.575068473815918}, 250.0: {'size': 57, 'filteredSize': 18, 'cost': 36.85774647887324, 'mem_MB': 1481.19921875, 'time_s': 4.840644598007202}, 125.0: {'size': 84, 'filteredSize': 23, 'cost': 41.909675713782164, 'mem_MB': 1484.828125, 'time_s': 4.962235689163208}}}, 0.3: {2: {4000: {'size': 1, 'filteredSize': 0, 'cost': 128.0, 'mem_MB': 19.64453125, 'time_s': 0.07669973373413086}, 2000.0: {'size': 2, 'filteredSize': 0, 'cost': 128.0, 'mem_MB': 233.23046875, 'time_s': 0.540287971496582}, 1000.0: {'size': 7, 'filteredSize': 0, 'cost': 128.0, 'mem_MB': 421.25, 'time_s': 0.9718663692474365}, 500.0: {'size': 22, 'filteredSize': 0, 'cost': 128.0, 'mem_MB': 894.00390625, 'time_s': 2.6226186752319336}, 250.0: {'size': 63, 'filteredSize': 0, 'cost': 128.0, 'mem_MB': 1448.51171875, 'time_s': 4.685434103012085}, 125.0: {'size': 99, 'filteredSize': 0, 'cost': 128.0, 'mem_MB': 1447.61328125, 'time_s': 4.767988681793213}}, 4: {4000: {'size': 1, 'filteredSize': 0, 'cost': 128.0, 'mem_MB': 19.640625, 'time_s': 0.07848143577575684}, 2000.0: {'size': 2, 'filteredSize': 0, 'cost': 128.0, 'mem_MB': 233.21875, 'time_s': 0.5306906700134277}, 1000.0: {'size': 7, 'filteredSize': 0, 'cost': 128.0, 'mem_MB': 422.0, 'time_s': 0.9626336097717285}, 500.0: {'size': 23, 'filteredSize': 0, 'cost': 128.0, 'mem_MB': 908.1015625, 'time_s': 2.6523780822753906}, 250.0: {'size': 55, 'filteredSize': 0, 'cost': 128.0, 'mem_MB': 1469.7265625, 'time_s': 4.704073905944824}, 125.0: {'size': 89, 'filteredSize': 0, 'cost': 128.0, 'mem_MB': 1470.859375, 'time_s': 4.870046377182007}}, 8: {4000: {'size': 1, 'filteredSize': 0, 'cost': 128.0, 'mem_MB': 19.64453125, 'time_s': 0.07740640640258789}, 2000.0: {'size': 2, 'filteredSize': 0, 'cost': 128.0, 'mem_MB': 233.2265625, 'time_s': 0.5315263271331787}, 1000.0: {'size': 7, 'filteredSize': 0, 'cost': 128.0, 'mem_MB': 421.57421875, 'time_s': 0.9598772525787354}, 500.0: {'size': 25, 'filteredSize': 0, 'cost': 128.0, 'mem_MB': 910.44921875, 'time_s': 2.555253267288208}, 250.0: {'size': 54, 'filteredSize': 0, 'cost': 128.0, 'mem_MB': 1481.75390625, 'time_s': 4.6191582679748535}, 125.0: {'size': 87, 'filteredSize': 0, 'cost': 128.0, 'mem_MB': 1483.17578125, 'time_s': 4.8185715675354}}, 16: {4000: {'size': 1, 'filteredSize': 0, 'cost': 128.0, 'mem_MB': 19.64453125, 'time_s': 0.07784056663513184}, 2000.0: {'size': 2, 'filteredSize': 0, 'cost': 128.0, 'mem_MB': 233.1953125, 'time_s': 0.5318508148193359}, 1000.0: {'size': 7, 'filteredSize': 0, 'cost': 128.0, 'mem_MB': 421.71484375, 'time_s': 0.9762487411499023}, 500.0: {'size': 25, 'filteredSize': 0, 'cost': 128.0, 'mem_MB': 910.953125, 'time_s': 2.5580482482910156}, 250.0: {'size': 57, 'filteredSize': 0, 'cost': 128.0, 'mem_MB': 1482.12109375, 'time_s': 4.803645133972168}, 125.0: {'size': 84, 'filteredSize': 0, 'cost': 128.0, 'mem_MB': 1485.3515625, 'time_s': 4.919614553451538}}, 32: {4000: {'size': 1, 'filteredSize': 0, 'cost': 128.0, 'mem_MB': 19.64453125, 'time_s': 0.07763314247131348}, 2000.0: {'size': 2, 'filteredSize': 0, 'cost': 128.0, 'mem_MB': 233.19921875, 'time_s': 0.5256521701812744}, 1000.0: {'size': 7, 'filteredSize': 0, 'cost': 128.0, 'mem_MB': 421.6953125, 'time_s': 0.974492073059082}, 500.0: {'size': 25, 'filteredSize': 0, 'cost': 128.0, 'mem_MB': 911.66796875, 'time_s': 2.5720267295837402}, 250.0: {'size': 57, 'filteredSize': 0, 'cost': 128.0, 'mem_MB': 1482.3046875, 'time_s': 4.861089706420898}, 125.0: {'size': 84, 'filteredSize': 0, 'cost': 128.0, 'mem_MB': 1485.5859375, 'time_s': 4.9071044921875}}, 64: {4000: {'size': 1, 'filteredSize': 0, 'cost': 128.0, 'mem_MB': 19.640625, 'time_s': 0.07813644409179688}, 2000.0: {'size': 2, 'filteredSize': 0, 'cost': 128.0, 'mem_MB': 233.19140625, 'time_s': 0.5250418186187744}, 1000.0: {'size': 7, 'filteredSize': 0, 'cost': 128.0, 'mem_MB': 422.69921875, 'time_s': 0.9712221622467041}, 500.0: {'size': 25, 'filteredSize': 0, 'cost': 128.0, 'mem_MB': 910.546875, 'time_s': 2.574735403060913}, 250.0: {'size': 57, 'filteredSize': 0, 'cost': 128.0, 'mem_MB': 1481.6875, 'time_s': 4.821501970291138}, 125.0: {'size': 84, 'filteredSize': 0, 'cost': 128.0, 'mem_MB': 1485.65625, 'time_s': 4.931427001953125}}}}
# ----------------------------
# Read LaTeX file
# ----------------------------
latex = latex_path.read_text()

# Extract all tabular blocks safely
tabular_blocks = re.findall(
    r'\\begin\{tabular\}.*?\\end\{tabular\}',
    latex,
    flags=re.DOTALL
)

# ----------------------------
# ✅ Top 2 per delta → Global Top 5
# ----------------------------
def select_top2_per_delta_then_top5(c2_data):
    delta_groups = {}

    for l, delta_dict in c2_data.items():
        for delta, info in delta_dict.items():
            delta_groups.setdefault(delta, []).append((l, delta, info))

    selected = []
    for delta, entries in delta_groups.items():
        entries.sort(key=lambda x: x[2]["cost"])
        selected.extend(entries[:1])

    selected.sort(key=lambda x: x[2]["cost"])
    return selected[:5]

# ----------------------------
# Format rows
# ----------------------------
def make_rows(entries, c2_vec):
    rows = []

    for l, delta, info in entries:
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

    # ✅ Always add 2 placeholder rows
    #placeholder = (
    #    f"basic-na & -inf & nan & {c2_vec} & "
    #    f"0 & -inf & nan & 128.00 & -inf & nan\\\\"
    #)
    #rows.append(placeholder)
    #rows.append(placeholder)

    return rows

# ----------------------------
# Update each tabular WITHOUT dropping old rows
# ----------------------------
updated_blocks = []

for block in tabular_blocks:

    match = re.search(r'\((?:1\.0,\s*([0-9.]+),\s*128\.0)', block)
    if not match:
        updated_blocks.append(block)
        continue

    c2_val = float(match.group(1))
    if c2_val not in data:
        updated_blocks.append(block)
        continue

    c2_vec = f"(1.0, {c2_val}, 128.0)"

    top5 = select_top2_per_delta_then_top5(data[c2_val])
    new_rows = make_rows(top5, c2_vec)

    lines = block.splitlines()

    bottom_idx = next(i for i, l in enumerate(lines) if r"\bottomrule" in l)

    # ✅ Insert BEFORE bottomrule (keep original rows!)
    new_block = "\n".join(
        lines[:bottom_idx] +
        new_rows +
        lines[bottom_idx:]
    )

    updated_blocks.append(new_block)

# ----------------------------
# Replace in full LaTeX
# ----------------------------
for old, new in zip(tabular_blocks, updated_blocks):
    latex = latex.replace(old, new)

# ----------------------------
# Write output
# ----------------------------
output_path.write_text(latex)

print(f"\n✅ SAFE UPDATE COMPLETE")
print(f"✅ Original rows preserved")
print(f"✅ New rows appended")
print(f"✅ Output written to:")
print(output_path)
