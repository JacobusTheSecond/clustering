import matplotlib.pyplot as plt
import json
import os 
import numpy as np

file_path = "result.json"

result = {}

measure = "macroF1"

measures = {"acc" : "Accuracy",
            "macroPrec" : "Macro-average precision",
            "trueAcc" : "Accuracy without unlabeled segments",
            "macroF1" : "Macro-average F1",
            "macroRec" : "Macro-average recall",}


if os.path.exists(file_path):
    with open(file_path, "r") as file:
        result = json.load(file)

colors = {
    "KlCluster": 'blue',
    "gmm": 'orange',
    "aca": 'green',
    "haca": 'red',
    "tmm"  : 'pink'
}
x = list(range(int(min(result.keys())), max([int(key) for key in result.keys()])+1))

# Acc values
acc_values = {method: [result[str(i)][measure][method] for i in x] for method in colors.keys()}

# Define colors for each method

# Converting line plot to bar chart
plt.figure(figsize=(12, 7))
bar_width = 0.15  # Width of each bar
index = np.arange(len(x))

for i, (method, values) in enumerate(acc_values.items()):
    plt.bar(index + i * bar_width, values, bar_width, label=method, color=colors[method])

plt.xlabel('TAG')
plt.ylabel(measures[measure])
plt.title(f'{measures[measure]} Values for Different Methods')
plt.xticks(index + bar_width, labels=x)
plt.legend()
plt.grid(True)
plt.show()