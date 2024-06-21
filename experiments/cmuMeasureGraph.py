import matplotlib.pyplot as plt
import json
import os 

file_path = "result.json"

result = {}

measure = "macroPrec"

measures = {"acc" : "Accuracy",
            "macroPrec" : "Macro-average precision",
            "macroRec" : "Macro-average recall",}


if os.path.exists(file_path):
    with open(file_path, "r") as file:
        result = json.load(file)

colors = {
    "KlCluster": 'blue',
    "gmm": 'orange',
    "aca": 'green',
    "haca": 'red'
}
x = list(range(int(min(result.keys())), max([int(key) for key in result.keys()])+1))

# Acc values
acc_values = {method: [result[str(i)][measure][method] for i in x] for method in colors.keys()}

# Define colors for each method


# Plotting
plt.figure(figsize=(10, 6))
for method, values in acc_values.items():
    plt.plot(x, values, label=method, marker='o', color=colors[method])

plt.xlabel('TAG')
plt.ylabel('Accuracy')
plt.title(f'{measures[measure]} Values for Different Methods')
plt.xticks(x, labels=x)
plt.legend()
plt.grid(True)
plt.show()