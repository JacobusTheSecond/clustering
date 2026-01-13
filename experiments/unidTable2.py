import re
from pathlib import Path


# ----------------------------
# ✅ Data format:
# data[c2][l][delta] = metrics
# ----------------------------
data = {3: {2: {32: {'size': 3, 'filteredSize': 3, 'cost': 99.0, 'mem_MB': 3188.75390625, 'time_s': 4.756692171096802}, 16.0: {'size': 9, 'filteredSize': 8, 'cost': 56.1724808012108, 'mem_MB': 4034.49609375, 'time_s': 7.951675891876221}, 8.0: {'size': 31, 'filteredSize': 19, 'cost': 47.48957379622191, 'mem_MB': 2446.33203125, 'time_s': 7.107620000839233}, 4.0: {'size': 105, 'filteredSize': 53, 'cost': 78.98954576879238, 'mem_MB': 1209.59375, 'time_s': 5.070381164550781}, 2.0: {'size': 246, 'filteredSize': 75, 'cost': 130.52228180646125, 'mem_MB': 702.0, 'time_s': 4.186448574066162}, 1.0: {'size': 455, 'filteredSize': 94, 'cost': 188.55338897413827, 'mem_MB': 448.48046875, 'time_s': 4.888254880905151}}, 4: {32: {'size': 3, 'filteredSize': 3, 'cost': 99.0, 'mem_MB': 3187.86328125, 'time_s': 4.744678974151611}, 16.0: {'size': 9, 'filteredSize': 8, 'cost': 56.1724808012108, 'mem_MB': 4035.8046875, 'time_s': 7.940945386886597}, 8.0: {'size': 31, 'filteredSize': 19, 'cost': 47.48957379622191, 'mem_MB': 2446.89453125, 'time_s': 7.149437427520752}, 4.0: {'size': 105, 'filteredSize': 52, 'cost': 78.75387245651075, 'mem_MB': 1209.703125, 'time_s': 5.031826972961426}, 2.0: {'size': 249, 'filteredSize': 79, 'cost': 129.1922868513986, 'mem_MB': 703.46484375, 'time_s': 4.1434006690979}, 1.0: {'size': 423, 'filteredSize': 96, 'cost': 182.54119727534714, 'mem_MB': 449.8359375, 'time_s': 4.729304313659668}}, 8: {32: {'size': 3, 'filteredSize': 3, 'cost': 99.0, 'mem_MB': 3188.875, 'time_s': 4.747220277786255}, 16.0: {'size': 9, 'filteredSize': 8, 'cost': 56.1724808012108, 'mem_MB': 4034.296875, 'time_s': 8.013218641281128}, 8.0: {'size': 31, 'filteredSize': 19, 'cost': 47.48957379622191, 'mem_MB': 2448.00390625, 'time_s': 7.166399955749512}, 4.0: {'size': 105, 'filteredSize': 52, 'cost': 78.75387245651075, 'mem_MB': 1209.59375, 'time_s': 5.035888433456421}, 2.0: {'size': 249, 'filteredSize': 79, 'cost': 129.1922868513986, 'mem_MB': 702.37109375, 'time_s': 4.278985023498535}, 1.0: {'size': 414, 'filteredSize': 97, 'cost': 181.22103184999438, 'mem_MB': 449.6796875, 'time_s': 4.581095933914185}}, 16: {32: {'size': 3, 'filteredSize': 3, 'cost': 99.0, 'mem_MB': 3187.19921875, 'time_s': 4.751485824584961}, 16.0: {'size': 9, 'filteredSize': 8, 'cost': 56.1724808012108, 'mem_MB': 4036.1953125, 'time_s': 8.00793170928955}, 8.0: {'size': 31, 'filteredSize': 19, 'cost': 47.48957379622191, 'mem_MB': 2446.76953125, 'time_s': 7.069660902023315}, 4.0: {'size': 105, 'filteredSize': 52, 'cost': 78.75387245651075, 'mem_MB': 1209.59765625, 'time_s': 5.125584602355957}, 2.0: {'size': 249, 'filteredSize': 79, 'cost': 129.1922868513986, 'mem_MB': 702.31640625, 'time_s': 4.259223699569702}, 1.0: {'size': 414, 'filteredSize': 97, 'cost': 181.22103184999438, 'mem_MB': 450.875, 'time_s': 4.695981025695801}}, 32: {32: {'size': 3, 'filteredSize': 3, 'cost': 99.0, 'mem_MB': 3190.3671875, 'time_s': 4.737264156341553}, 16.0: {'size': 9, 'filteredSize': 8, 'cost': 56.1724808012108, 'mem_MB': 4035.31640625, 'time_s': 8.030378103256226}, 8.0: {'size': 31, 'filteredSize': 19, 'cost': 47.48957379622191, 'mem_MB': 2448.39453125, 'time_s': 7.155622720718384}, 4.0: {'size': 105, 'filteredSize': 52, 'cost': 78.75387245651075, 'mem_MB': 1209.14453125, 'time_s': 5.235355854034424}, 2.0: {'size': 249, 'filteredSize': 79, 'cost': 129.1922868513986, 'mem_MB': 701.73046875, 'time_s': 4.40920877456665}, 1.0: {'size': 414, 'filteredSize': 97, 'cost': 181.22103184999438, 'mem_MB': 448.62890625, 'time_s': 4.616863012313843}}, 64: {32: {'size': 3, 'filteredSize': 3, 'cost': 99.0, 'mem_MB': 3187.8203125, 'time_s': 4.774210214614868}, 16.0: {'size': 9, 'filteredSize': 8, 'cost': 56.1724808012108, 'mem_MB': 4034.2890625, 'time_s': 7.947916030883789}, 8.0: {'size': 31, 'filteredSize': 19, 'cost': 47.48957379622191, 'mem_MB': 2447.078125, 'time_s': 7.188919544219971}, 4.0: {'size': 105, 'filteredSize': 52, 'cost': 78.75387245651075, 'mem_MB': 1209.3359375, 'time_s': 5.084000587463379}, 2.0: {'size': 249, 'filteredSize': 79, 'cost': 129.1922868513986, 'mem_MB': 700.9375, 'time_s': 4.4281415939331055}, 1.0: {'size': 414, 'filteredSize': 97, 'cost': 181.22103184999438, 'mem_MB': 449.15625, 'time_s': 4.549536228179932}}}, 30: {2: {32: {'size': 3, 'filteredSize': 0, 'cost': 362.0, 'mem_MB': 3188.41015625, 'time_s': 4.749939918518066}, 16.0: {'size': 9, 'filteredSize': 0, 'cost': 362.0, 'mem_MB': 4035.43359375, 'time_s': 8.063246965408325}, 8.0: {'size': 31, 'filteredSize': 0, 'cost': 362.0, 'mem_MB': 2445.953125, 'time_s': 7.020982027053833}, 4.0: {'size': 105, 'filteredSize': 0, 'cost': 362.0, 'mem_MB': 1209.03515625, 'time_s': 4.8359410762786865}, 2.0: {'size': 246, 'filteredSize': 0, 'cost': 362.0, 'mem_MB': 701.32421875, 'time_s': 3.0127487182617188}, 1.0: {'size': 455, 'filteredSize': 0, 'cost': 362.0, 'mem_MB': 446.890625, 'time_s': 2.2558159828186035}}, 4: {32: {'size': 3, 'filteredSize': 0, 'cost': 362.0, 'mem_MB': 3187.7421875, 'time_s': 4.761719465255737}, 16.0: {'size': 9, 'filteredSize': 0, 'cost': 362.0, 'mem_MB': 4035.546875, 'time_s': 7.924484491348267}, 8.0: {'size': 31, 'filteredSize': 0, 'cost': 362.0, 'mem_MB': 2447.0703125, 'time_s': 7.108106851577759}, 4.0: {'size': 105, 'filteredSize': 0, 'cost': 362.0, 'mem_MB': 1209.109375, 'time_s': 4.6970601081848145}, 2.0: {'size': 249, 'filteredSize': 0, 'cost': 362.0, 'mem_MB': 701.28515625, 'time_s': 3.095047950744629}, 1.0: {'size': 423, 'filteredSize': 0, 'cost': 362.0, 'mem_MB': 450.53125, 'time_s': 2.3565168380737305}}, 8: {32: {'size': 3, 'filteredSize': 0, 'cost': 362.0, 'mem_MB': 3187.3828125, 'time_s': 4.7244651317596436}, 16.0: {'size': 9, 'filteredSize': 0, 'cost': 362.0, 'mem_MB': 4035.515625, 'time_s': 8.07278060913086}, 8.0: {'size': 31, 'filteredSize': 0, 'cost': 362.0, 'mem_MB': 2447.359375, 'time_s': 7.170115947723389}, 4.0: {'size': 105, 'filteredSize': 0, 'cost': 362.0, 'mem_MB': 1209.328125, 'time_s': 4.6834399700164795}, 2.0: {'size': 249, 'filteredSize': 0, 'cost': 362.0, 'mem_MB': 702.1640625, 'time_s': 3.0750844478607178}, 1.0: {'size': 414, 'filteredSize': 0, 'cost': 362.0, 'mem_MB': 448.84765625, 'time_s': 2.3516592979431152}}, 16: {32: {'size': 3, 'filteredSize': 0, 'cost': 362.0, 'mem_MB': 3187.8671875, 'time_s': 4.745477199554443}, 16.0: {'size': 9, 'filteredSize': 0, 'cost': 362.0, 'mem_MB': 4036.4453125, 'time_s': 7.970113515853882}, 8.0: {'size': 31, 'filteredSize': 0, 'cost': 362.0, 'mem_MB': 2447.484375, 'time_s': 7.101258039474487}, 4.0: {'size': 105, 'filteredSize': 0, 'cost': 362.0, 'mem_MB': 1208.9296875, 'time_s': 4.564249515533447}, 2.0: {'size': 249, 'filteredSize': 0, 'cost': 362.0, 'mem_MB': 702.0078125, 'time_s': 2.893191337585449}, 1.0: {'size': 414, 'filteredSize': 0, 'cost': 362.0, 'mem_MB': 449.8359375, 'time_s': 2.2220299243927}}, 32: {32: {'size': 3, 'filteredSize': 0, 'cost': 362.0, 'mem_MB': 3187.3203125, 'time_s': 4.764532566070557}, 16.0: {'size': 9, 'filteredSize': 0, 'cost': 362.0, 'mem_MB': 4035.125, 'time_s': 7.998162746429443}, 8.0: {'size': 31, 'filteredSize': 0, 'cost': 362.0, 'mem_MB': 2446.03515625, 'time_s': 7.080420970916748}, 4.0: {'size': 105, 'filteredSize': 0, 'cost': 362.0, 'mem_MB': 1209.609375, 'time_s': 4.761027574539185}, 2.0: {'size': 249, 'filteredSize': 0, 'cost': 362.0, 'mem_MB': 703.46875, 'time_s': 2.9533936977386475}, 1.0: {'size': 414, 'filteredSize': 0, 'cost': 362.0, 'mem_MB': 449.5859375, 'time_s': 2.426856756210327}}, 64: {32: {'size': 3, 'filteredSize': 0, 'cost': 362.0, 'mem_MB': 3188.19921875, 'time_s': 4.759751081466675}, 16.0: {'size': 9, 'filteredSize': 0, 'cost': 362.0, 'mem_MB': 4034.76171875, 'time_s': 7.895925283432007}, 8.0: {'size': 31, 'filteredSize': 0, 'cost': 362.0, 'mem_MB': 2446.5, 'time_s': 7.0188517570495605}, 4.0: {'size': 105, 'filteredSize': 0, 'cost': 362.0, 'mem_MB': 1209.57421875, 'time_s': 4.70081901550293}, 2.0: {'size': 249, 'filteredSize': 0, 'cost': 362.0, 'mem_MB': 701.17578125, 'time_s': 3.02068829536438}, 1.0: {'size': 414, 'filteredSize': 0, 'cost': 362.0, 'mem_MB': 450.4375, 'time_s': 2.177095413208008}}}, 300: {2: {32: {'size': 3, 'filteredSize': 0, 'cost': 362.0, 'mem_MB': 3188.140625, 'time_s': 4.708207130432129}, 16.0: {'size': 9, 'filteredSize': 0, 'cost': 362.0, 'mem_MB': 4035.40625, 'time_s': 7.977485418319702}, 8.0: {'size': 31, 'filteredSize': 0, 'cost': 362.0, 'mem_MB': 2446.8203125, 'time_s': 7.146981477737427}, 4.0: {'size': 105, 'filteredSize': 0, 'cost': 362.0, 'mem_MB': 1209.58203125, 'time_s': 4.682678937911987}, 2.0: {'size': 246, 'filteredSize': 0, 'cost': 362.0, 'mem_MB': 701.2578125, 'time_s': 3.0172441005706787}, 1.0: {'size': 455, 'filteredSize': 0, 'cost': 362.0, 'mem_MB': 447.97265625, 'time_s': 2.392378330230713}}, 4: {32: {'size': 3, 'filteredSize': 0, 'cost': 362.0, 'mem_MB': 3187.97265625, 'time_s': 4.775926828384399}, 16.0: {'size': 9, 'filteredSize': 0, 'cost': 362.0, 'mem_MB': 4036.6171875, 'time_s': 7.982968091964722}, 8.0: {'size': 31, 'filteredSize': 0, 'cost': 362.0, 'mem_MB': 2446.54296875, 'time_s': 7.071635961532593}, 4.0: {'size': 105, 'filteredSize': 0, 'cost': 362.0, 'mem_MB': 1209.3515625, 'time_s': 4.686020851135254}, 2.0: {'size': 249, 'filteredSize': 0, 'cost': 362.0, 'mem_MB': 702.65625, 'time_s': 3.1504926681518555}, 1.0: {'size': 423, 'filteredSize': 0, 'cost': 362.0, 'mem_MB': 449.125, 'time_s': 2.270979166030884}}, 8: {32: {'size': 3, 'filteredSize': 0, 'cost': 362.0, 'mem_MB': 3188.25390625, 'time_s': 4.7113707065582275}, 16.0: {'size': 9, 'filteredSize': 0, 'cost': 362.0, 'mem_MB': 4034.30078125, 'time_s': 7.900886535644531}, 8.0: {'size': 31, 'filteredSize': 0, 'cost': 362.0, 'mem_MB': 2446.8125, 'time_s': 7.066399097442627}, 4.0: {'size': 105, 'filteredSize': 0, 'cost': 362.0, 'mem_MB': 1209.15625, 'time_s': 4.609228849411011}, 2.0: {'size': 249, 'filteredSize': 0, 'cost': 362.0, 'mem_MB': 700.8515625, 'time_s': 3.082979440689087}, 1.0: {'size': 414, 'filteredSize': 0, 'cost': 362.0, 'mem_MB': 449.375, 'time_s': 2.205721855163574}}, 16: {32: {'size': 3, 'filteredSize': 0, 'cost': 362.0, 'mem_MB': 3187.9296875, 'time_s': 4.7576377391815186}, 16.0: {'size': 9, 'filteredSize': 0, 'cost': 362.0, 'mem_MB': 4035.0234375, 'time_s': 7.922630310058594}, 8.0: {'size': 31, 'filteredSize': 0, 'cost': 362.0, 'mem_MB': 2445.8671875, 'time_s': 7.025359392166138}, 4.0: {'size': 105, 'filteredSize': 0, 'cost': 362.0, 'mem_MB': 1209.2109375, 'time_s': 4.64506459236145}, 2.0: {'size': 249, 'filteredSize': 0, 'cost': 362.0, 'mem_MB': 701.81640625, 'time_s': 3.075343132019043}, 1.0: {'size': 414, 'filteredSize': 0, 'cost': 362.0, 'mem_MB': 449.7265625, 'time_s': 2.3439836502075195}}, 32: {32: {'size': 3, 'filteredSize': 0, 'cost': 362.0, 'mem_MB': 3187.1796875, 'time_s': 4.721126556396484}, 16.0: {'size': 9, 'filteredSize': 0, 'cost': 362.0, 'mem_MB': 4034.52734375, 'time_s': 7.98747706413269}, 8.0: {'size': 31, 'filteredSize': 0, 'cost': 362.0, 'mem_MB': 2446.66015625, 'time_s': 7.138868093490601}, 4.0: {'size': 105, 'filteredSize': 0, 'cost': 362.0, 'mem_MB': 1209.53125, 'time_s': 4.847756862640381}, 2.0: {'size': 249, 'filteredSize': 0, 'cost': 362.0, 'mem_MB': 702.39453125, 'time_s': 3.0513248443603516}, 1.0: {'size': 414, 'filteredSize': 0, 'cost': 362.0, 'mem_MB': 450.7890625, 'time_s': 2.234740972518921}}, 64: {32: {'size': 3, 'filteredSize': 0, 'cost': 362.0, 'mem_MB': 3187.6015625, 'time_s': 4.716050386428833}, 16.0: {'size': 9, 'filteredSize': 0, 'cost': 362.0, 'mem_MB': 4035.29296875, 'time_s': 7.947446823120117}, 8.0: {'size': 31, 'filteredSize': 0, 'cost': 362.0, 'mem_MB': 2446.8828125, 'time_s': 7.21895956993103}, 4.0: {'size': 105, 'filteredSize': 0, 'cost': 362.0, 'mem_MB': 1209.53125, 'time_s': 4.708579778671265}, 2.0: {'size': 249, 'filteredSize': 0, 'cost': 362.0, 'mem_MB': 700.859375, 'time_s': 3.0887985229492188}, 1.0: {'size': 414, 'filteredSize': 0, 'cost': 362.0, 'mem_MB': 448.58984375, 'time_s': 2.3575994968414307}}}}
# ----------------------------
# Read LaTeX file
# ----------------------------

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

for c2_val in [3,30,300]:

    c2_vec = f"(1.0, {c2_val}, 362.0)"

    top5 = select_top2_per_delta_then_top5(data[c2_val])
    new_rows = make_rows(top5, c2_vec)
    
    print("\n".join(new_rows))
    print("-----")
# ----------------------------

