import re
from pathlib import Path

# Base folder of the script
script_dir = Path(__file__).parent

# Path to the LaTeX file relative to the script
latex_path = script_dir / "../discrete-subtrajectory-clustering-test-suite/data/plots/socg_table_athens_centers.tex"
latex_path = latex_path.resolve()
output_path = latex_path.with_name("socg_table_athens_centers_updated.tex")

# Your data: c2 -> l -> delta -> metrics
data = {0.003: {2: {4000: {'size': 1, 'filteredSize': 1, 'cost': 13.0, 'mem_MB': 19.640625, 'time_s': 0.07750797271728516}, 2000.0: {'size': 2, 'filteredSize': 2, 'cost': 8.0, 'mem_MB': 233.2109375, 'time_s': 0.5229527950286865}, 1000.0: {'size': 7, 'filteredSize': 4, 'cost': 7.766197183098598, 'mem_MB': 421.53515625, 'time_s': 0.950340986251831}, 500.0: {'size': 25, 'filteredSize': 11, 'cost': 15.790140845070425, 'mem_MB': 910.8203125, 'time_s': 2.551877021789551}, 250.0: {'size': 57, 'filteredSize': 18, 'cost': 30.10774647887324, 'mem_MB': 1481.19921875, 'time_s': 4.697320938110352}, 125.0: {'size': 84, 'filteredSize': 23, 'cost': 38.534675713782164, 'mem_MB': 1484.171875, 'time_s': 4.9354681968688965}}, 4: {4000: {'size': 1, 'filteredSize': 1, 'cost': 13.0, 'mem_MB': 19.6484375, 'time_s': 0.0884547233581543}, 2000.0: {'size': 2, 'filteredSize': 2, 'cost': 8.0, 'mem_MB': 233.21875, 'time_s': 0.5229661464691162}, 1000.0: {'size': 7, 'filteredSize': 4, 'cost': 7.766197183098598, 'mem_MB': 422.234375, 'time_s': 0.9753215312957764}, 500.0: {'size': 25, 'filteredSize': 11, 'cost': 15.790140845070425, 'mem_MB': 911.046875, 'time_s': 2.5654923915863037}, 250.0: {'size': 57, 'filteredSize': 18, 'cost': 30.10774647887324, 'mem_MB': 1481.8046875, 'time_s': 4.695620775222778}, 125.0: {'size': 84, 'filteredSize': 23, 'cost': 38.534675713782164, 'mem_MB': 1484.83203125, 'time_s': 4.956231117248535}}, 8: {4000: {'size': 1, 'filteredSize': 1, 'cost': 13.0, 'mem_MB': 19.64453125, 'time_s': 0.07901191711425781}, 2000.0: {'size': 2, 'filteredSize': 2, 'cost': 8.0, 'mem_MB': 233.2421875, 'time_s': 0.5235781669616699}, 1000.0: {'size': 7, 'filteredSize': 4, 'cost': 7.766197183098598, 'mem_MB': 421.5859375, 'time_s': 0.953432559967041}, 500.0: {'size': 25, 'filteredSize': 11, 'cost': 15.790140845070425, 'mem_MB': 911.109375, 'time_s': 2.6361029148101807}, 250.0: {'size': 57, 'filteredSize': 18, 'cost': 30.10774647887324, 'mem_MB': 1481.87890625, 'time_s': 4.8348708152771}, 125.0: {'size': 84, 'filteredSize': 23, 'cost': 38.534675713782164, 'mem_MB': 1485.625, 'time_s': 5.015235424041748}}, 16: {4000: {'size': 1, 'filteredSize': 1, 'cost': 13.0, 'mem_MB': 19.62109375, 'time_s': 0.07621240615844727}, 2000.0: {'size': 2, 'filteredSize': 2, 'cost': 8.0, 'mem_MB': 233.2421875, 'time_s': 0.5256426334381104}, 1000.0: {'size': 7, 'filteredSize': 4, 'cost': 7.766197183098598, 'mem_MB': 422.30078125, 'time_s': 0.9608607292175293}, 500.0: {'size': 25, 'filteredSize': 11, 'cost': 15.790140845070425, 'mem_MB': 911.26171875, 'time_s': 2.5859649181365967}, 250.0: {'size': 57, 'filteredSize': 18, 'cost': 30.10774647887324, 'mem_MB': 1481.5234375, 'time_s': 4.813166618347168}, 125.0: {'size': 84, 'filteredSize': 23, 'cost': 38.534675713782164, 'mem_MB': 1486.06640625, 'time_s': 4.9662840366363525}}, 32: {4000: {'size': 1, 'filteredSize': 1, 'cost': 13.0, 'mem_MB': 19.62109375, 'time_s': 0.07810473442077637}, 2000.0: {'size': 2, 'filteredSize': 2, 'cost': 8.0, 'mem_MB': 233.20703125, 'time_s': 0.5308337211608887}, 1000.0: {'size': 7, 'filteredSize': 4, 'cost': 7.766197183098598, 'mem_MB': 421.8515625, 'time_s': 0.9689781665802002}, 500.0: {'size': 25, 'filteredSize': 11, 'cost': 15.790140845070425, 'mem_MB': 910.515625, 'time_s': 2.5909998416900635}, 250.0: {'size': 57, 'filteredSize': 18, 'cost': 30.10774647887324, 'mem_MB': 1480.6484375, 'time_s': 4.894773006439209}, 125.0: {'size': 84, 'filteredSize': 23, 'cost': 38.534675713782164, 'mem_MB': 1484.78515625, 'time_s': 4.921672344207764}}, 64: {4000: {'size': 1, 'filteredSize': 1, 'cost': 13.0, 'mem_MB': 19.640625, 'time_s': 0.07657146453857422}, 2000.0: {'size': 2, 'filteredSize': 2, 'cost': 8.0, 'mem_MB': 233.2421875, 'time_s': 0.523306131362915}, 1000.0: {'size': 7, 'filteredSize': 4, 'cost': 7.766197183098598, 'mem_MB': 422.375, 'time_s': 0.9496965408325195}, 500.0: {'size': 25, 'filteredSize': 11, 'cost': 15.790140845070425, 'mem_MB': 910.51953125, 'time_s': 2.6143839359283447}, 250.0: {'size': 57, 'filteredSize': 18, 'cost': 30.10774647887324, 'mem_MB': 1482.31640625, 'time_s': 4.76694393157959}, 125.0: {'size': 84, 'filteredSize': 23, 'cost': 38.534675713782164, 'mem_MB': 1484.828125, 'time_s': 4.887715101242065}}}, 0.03: {2: {4000: {'size': 1, 'filteredSize': 1, 'cost': 121.0, 'mem_MB': 19.64453125, 'time_s': 0.07229113578796387}, 2000.0: {'size': 2, 'filteredSize': 2, 'cost': 62.0, 'mem_MB': 233.2421875, 'time_s': 0.5138039588928223}, 1000.0: {'size': 7, 'filteredSize': 4, 'cost': 34.7661971830986, 'mem_MB': 422.05859375, 'time_s': 0.9383499622344971}, 500.0: {'size': 25, 'filteredSize': 11, 'cost': 29.290140845070425, 'mem_MB': 911.4765625, 'time_s': 2.589209794998169}, 250.0: {'size': 57, 'filteredSize': 18, 'cost': 36.85774647887324, 'mem_MB': 1481.3359375, 'time_s': 4.822797536849976}, 125.0: {'size': 84, 'filteredSize': 23, 'cost': 41.909675713782164, 'mem_MB': 1485.27734375, 'time_s': 4.929623365402222}}, 4: {4000: {'size': 1, 'filteredSize': 1, 'cost': 121.0, 'mem_MB': 19.65234375, 'time_s': 0.07693648338317871}, 2000.0: {'size': 2, 'filteredSize': 2, 'cost': 62.0, 'mem_MB': 233.24609375, 'time_s': 0.5306172370910645}, 1000.0: {'size': 7, 'filteredSize': 4, 'cost': 34.7661971830986, 'mem_MB': 422.6640625, 'time_s': 0.9716897010803223}, 500.0: {'size': 25, 'filteredSize': 11, 'cost': 29.290140845070425, 'mem_MB': 911.0234375, 'time_s': 2.5863265991210938}, 250.0: {'size': 57, 'filteredSize': 18, 'cost': 36.85774647887324, 'mem_MB': 1481.9921875, 'time_s': 4.885124444961548}, 125.0: {'size': 84, 'filteredSize': 23, 'cost': 41.909675713782164, 'mem_MB': 1486.8046875, 'time_s': 4.901963472366333}}, 8: {4000: {'size': 1, 'filteredSize': 1, 'cost': 121.0, 'mem_MB': 19.640625, 'time_s': 0.07686734199523926}, 2000.0: {'size': 2, 'filteredSize': 2, 'cost': 62.0, 'mem_MB': 233.2265625, 'time_s': 0.5204510688781738}, 1000.0: {'size': 7, 'filteredSize': 4, 'cost': 34.7661971830986, 'mem_MB': 421.91796875, 'time_s': 0.9705729484558105}, 500.0: {'size': 25, 'filteredSize': 11, 'cost': 29.290140845070425, 'mem_MB': 910.640625, 'time_s': 2.5756609439849854}, 250.0: {'size': 57, 'filteredSize': 18, 'cost': 36.85774647887324, 'mem_MB': 1481.375, 'time_s': 4.861271142959595}, 125.0: {'size': 84, 'filteredSize': 23, 'cost': 41.909675713782164, 'mem_MB': 1485.2734375, 'time_s': 4.916283369064331}}, 16: {4000: {'size': 1, 'filteredSize': 1, 'cost': 121.0, 'mem_MB': 19.640625, 'time_s': 0.07815694808959961}, 2000.0: {'size': 2, 'filteredSize': 2, 'cost': 62.0, 'mem_MB': 233.24609375, 'time_s': 0.5160574913024902}, 1000.0: {'size': 7, 'filteredSize': 4, 'cost': 34.7661971830986, 'mem_MB': 421.5859375, 'time_s': 0.9341490268707275}, 500.0: {'size': 25, 'filteredSize': 11, 'cost': 29.290140845070425, 'mem_MB': 911.703125, 'time_s': 2.5818874835968018}, 250.0: {'size': 57, 'filteredSize': 18, 'cost': 36.85774647887324, 'mem_MB': 1480.7421875, 'time_s': 4.857981443405151}, 125.0: {'size': 84, 'filteredSize': 23, 'cost': 41.909675713782164, 'mem_MB': 1486.83984375, 'time_s': 4.859561204910278}}, 32: {4000: {'size': 1, 'filteredSize': 1, 'cost': 121.0, 'mem_MB': 19.640625, 'time_s': 0.07699084281921387}, 2000.0: {'size': 2, 'filteredSize': 2, 'cost': 62.0, 'mem_MB': 233.2421875, 'time_s': 0.5301308631896973}, 1000.0: {'size': 7, 'filteredSize': 4, 'cost': 34.7661971830986, 'mem_MB': 422.23828125, 'time_s': 0.9695925712585449}, 500.0: {'size': 25, 'filteredSize': 11, 'cost': 29.290140845070425, 'mem_MB': 911.0234375, 'time_s': 2.606442928314209}, 250.0: {'size': 57, 'filteredSize': 18, 'cost': 36.85774647887324, 'mem_MB': 1480.8125, 'time_s': 4.825195074081421}, 125.0: {'size': 84, 'filteredSize': 23, 'cost': 41.909675713782164, 'mem_MB': 1485.5546875, 'time_s': 4.894541501998901}}, 64: {4000: {'size': 1, 'filteredSize': 1, 'cost': 121.0, 'mem_MB': 19.6171875, 'time_s': 0.07768654823303223}, 2000.0: {'size': 2, 'filteredSize': 2, 'cost': 62.0, 'mem_MB': 233.2421875, 'time_s': 0.5223870277404785}, 1000.0: {'size': 7, 'filteredSize': 4, 'cost': 34.7661971830986, 'mem_MB': 422.8515625, 'time_s': 0.9424035549163818}, 500.0: {'size': 25, 'filteredSize': 11, 'cost': 29.290140845070425, 'mem_MB': 910.48828125, 'time_s': 2.5387959480285645}, 250.0: {'size': 57, 'filteredSize': 18, 'cost': 36.85774647887324, 'mem_MB': 1481.6796875, 'time_s': 4.799349546432495}, 125.0: {'size': 84, 'filteredSize': 23, 'cost': 41.909675713782164, 'mem_MB': 1484.2421875, 'time_s': 4.917693376541138}}}, 0.3: {2: {4000: {'size': 1, 'filteredSize': 0, 'cost': 128.0, 'mem_MB': 19.63671875, 'time_s': 0.07898068428039551}, 2000.0: {'size': 2, 'filteredSize': 0, 'cost': 128.0, 'mem_MB': 233.203125, 'time_s': 0.5281362533569336}, 1000.0: {'size': 7, 'filteredSize': 0, 'cost': 128.0, 'mem_MB': 421.5390625, 'time_s': 0.9482419490814209}, 500.0: {'size': 25, 'filteredSize': 0, 'cost': 128.0, 'mem_MB': 910.98828125, 'time_s': 2.617880344390869}, 250.0: {'size': 57, 'filteredSize': 0, 'cost': 128.0, 'mem_MB': 1481.421875, 'time_s': 4.826529026031494}, 125.0: {'size': 84, 'filteredSize': 0, 'cost': 128.0, 'mem_MB': 1484.21875, 'time_s': 4.889861822128296}}, 4: {4000: {'size': 1, 'filteredSize': 0, 'cost': 128.0, 'mem_MB': 19.640625, 'time_s': 0.07265758514404297}, 2000.0: {'size': 2, 'filteredSize': 0, 'cost': 128.0, 'mem_MB': 233.234375, 'time_s': 0.5290994644165039}, 1000.0: {'size': 7, 'filteredSize': 0, 'cost': 128.0, 'mem_MB': 422.46875, 'time_s': 0.9421172142028809}, 500.0: {'size': 25, 'filteredSize': 0, 'cost': 128.0, 'mem_MB': 911.07421875, 'time_s': 2.5762717723846436}, 250.0: {'size': 57, 'filteredSize': 0, 'cost': 128.0, 'mem_MB': 1481.68359375, 'time_s': 4.808232545852661}, 125.0: {'size': 84, 'filteredSize': 0, 'cost': 128.0, 'mem_MB': 1484.20703125, 'time_s': 4.8996992111206055}}, 8: {4000: {'size': 1, 'filteredSize': 0, 'cost': 128.0, 'mem_MB': 19.63671875, 'time_s': 0.07890033721923828}, 2000.0: {'size': 2, 'filteredSize': 0, 'cost': 128.0, 'mem_MB': 233.22265625, 'time_s': 0.5284388065338135}, 1000.0: {'size': 7, 'filteredSize': 0, 'cost': 128.0, 'mem_MB': 421.56640625, 'time_s': 0.9557778835296631}, 500.0: {'size': 25, 'filteredSize': 0, 'cost': 128.0, 'mem_MB': 911.3046875, 'time_s': 2.6221790313720703}, 250.0: {'size': 57, 'filteredSize': 0, 'cost': 128.0, 'mem_MB': 1481.9453125, 'time_s': 4.976045370101929}, 125.0: {'size': 84, 'filteredSize': 0, 'cost': 128.0, 'mem_MB': 1484.296875, 'time_s': 4.973398447036743}}, 16: {4000: {'size': 1, 'filteredSize': 0, 'cost': 128.0, 'mem_MB': 19.640625, 'time_s': 0.08063721656799316}, 2000.0: {'size': 2, 'filteredSize': 0, 'cost': 128.0, 'mem_MB': 233.21484375, 'time_s': 0.5289211273193359}, 1000.0: {'size': 7, 'filteredSize': 0, 'cost': 128.0, 'mem_MB': 422.7578125, 'time_s': 0.9623382091522217}, 500.0: {'size': 25, 'filteredSize': 0, 'cost': 128.0, 'mem_MB': 910.48828125, 'time_s': 2.594052314758301}, 250.0: {'size': 57, 'filteredSize': 0, 'cost': 128.0, 'mem_MB': 1480.59375, 'time_s': 4.728553295135498}, 125.0: {'size': 84, 'filteredSize': 0, 'cost': 128.0, 'mem_MB': 1485.234375, 'time_s': 4.972213268280029}}, 32: {4000: {'size': 1, 'filteredSize': 0, 'cost': 128.0, 'mem_MB': 19.6328125, 'time_s': 0.07616567611694336}, 2000.0: {'size': 2, 'filteredSize': 0, 'cost': 128.0, 'mem_MB': 233.22265625, 'time_s': 0.5323226451873779}, 1000.0: {'size': 7, 'filteredSize': 0, 'cost': 128.0, 'mem_MB': 421.56640625, 'time_s': 0.9630138874053955}, 500.0: {'size': 25, 'filteredSize': 0, 'cost': 128.0, 'mem_MB': 910.4765625, 'time_s': 2.5946362018585205}, 250.0: {'size': 57, 'filteredSize': 0, 'cost': 128.0, 'mem_MB': 1480.8046875, 'time_s': 4.91428279876709}, 125.0: {'size': 84, 'filteredSize': 0, 'cost': 128.0, 'mem_MB': 1485.23828125, 'time_s': 4.982137680053711}}, 64: {4000: {'size': 1, 'filteredSize': 0, 'cost': 128.0, 'mem_MB': 19.6171875, 'time_s': 0.07780575752258301}, 2000.0: {'size': 2, 'filteredSize': 0, 'cost': 128.0, 'mem_MB': 233.234375, 'time_s': 0.526909589767456}, 1000.0: {'size': 7, 'filteredSize': 0, 'cost': 128.0, 'mem_MB': 421.78515625, 'time_s': 0.9512476921081543}, 500.0: {'size': 25, 'filteredSize': 0, 'cost': 128.0, 'mem_MB': 910.7265625, 'time_s': 2.5985772609710693}, 250.0: {'size': 57, 'filteredSize': 0, 'cost': 128.0, 'mem_MB': 1481.66796875, 'time_s': 4.934208393096924}, 125.0: {'size': 84, 'filteredSize': 0, 'cost': 128.0, 'mem_MB': 1485.72265625, 'time_s': 4.953621864318848}}}}

# Read LaTeX file
with open(latex_path, "r") as f:
    latex = f.read()

# Split into individual tabulars
tabulars = re.split(r'(\\begin\{tabular\}.*?\\end\{tabular\})', latex, flags=re.DOTALL)
tabulars = [t for t in tabulars if t.strip().startswith(r"\begin{tabular}")]

# ----------------------------
# Helper: Select best per delta, then top 5
# ----------------------------
def select_unique_delta_top5(c2_data):
    best_per_delta = {}

    for l, deltas in c2_data.items():
        for delta, info in deltas.items():
            if delta not in best_per_delta or info["cost"] < best_per_delta[delta][2]["cost"]:
                best_per_delta[delta] = (l, delta, info)

    # Sort by cost and take top 5
    result = sorted(best_per_delta.values(), key=lambda x: x[2]["cost"])[:5]
    return result

# ----------------------------
# Helper: Generate LaTeX rows
# ----------------------------
def generate_rows(top5, c2_vec):
    rows = []

    for l, delta, info in top5:
        name = f"basic-{l}"
        seconds = info["time_s"]
        gbytes = info["mem_MB"] / 1024
        clusters = info["filteredSize"]
        max_frechet = delta
        avg_frechet = delta
        kcenters = info["cost"]
        max_size = ""#info["size"]
        avg_size = ""#info["size"]

        row = (
            f"{name} & {seconds:.2f} & {gbytes:.2f} & {c2_vec} & "
            f"{clusters} & {max_frechet} & {avg_frechet} & {kcenters:.2f} & "
            f"{max_size} & {avg_size}\\\\"
        )
        rows.append(row)

    # ✅ Two placeholder rows
    #placeholder = f"basic-na & -inf & nan & {c2_vec} & 0 & -inf & nan & 128.00 & -inf & nan\\\\"
    #rows.append(placeholder)
    #rows.append(placeholder)

    return rows

# ----------------------------
# Update each tabular
# ----------------------------
updated_tabulars = []

for tab in tabulars:
    # Detect c2
    match = re.search(r'\((?:1\.0,\s*([0-9.]+),\s*128\.0)', tab)
    if not match:
        updated_tabulars.append(tab)
        continue

    c2_val = float(match.group(1))
    if c2_val not in data:
        updated_tabulars.append(tab)
        continue

    c2_vec = (1.0, c2_val, 128.0)

    # ✅ Unique-delta best-cost selection
    top5 = select_unique_delta_top5(data[c2_val])

    lines = tab.splitlines()
    bottom_idx = next(i for i, line in enumerate(lines) if r"\bottomrule" in line)
    header_idx = next(i for i, line in enumerate(lines) if r"\midrule" in line)

    new_content = generate_rows(top5, c2_vec)

    new_tab = "\n".join(lines[:header_idx + 1] + new_content + lines[bottom_idx:])
    updated_tabulars.append(new_tab)

# ----------------------------
# Replace tabulars in LaTeX
# ----------------------------
for old, new in zip(tabulars, updated_tabulars):
    latex = latex.replace(old, new)

# Save output
with open(output_path, "w") as f:
    f.write(latex)

print(f"✅ Updated LaTeX table written to:\n{output_path}")
