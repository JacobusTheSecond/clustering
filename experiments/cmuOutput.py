import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import gmean


# Load the CSV file
file_path = 'output.csv'
df = pd.read_csv(file_path)

measure = 'Accuracy'

# Filter only the necessary columns: TAG, Method, 'measure', Complexity
df_filtered = df[['TAG', 'Method', measure, 'Complexity']]

# Separate the methods without complexity (Tmm, Haca, Aca, Gmm)
df_tmm = df_filtered[df_filtered['Method'] == 'Tmm']
df_aca = df_filtered[df_filtered['Method'] == 'Aca']
df_haca = df_filtered[df_filtered['Method'] == 'Haca']
df_gmm = df_filtered[df_filtered['Method'] == 'Gmm']

# For KlCluster, filter only the rows where complexity is 5, 10 or 15
df_klcluster1 = df_filtered[(df_filtered['Method'] == 'KlCluster') & (df_filtered['Complexity'].isin([5]))]
df_klcluster1['Method'].replace('KlCluster', 'KlCluster 5', inplace=True)

df_klcluster2 = df_filtered[(df_filtered['Method'] == 'KlCluster') & (df_filtered['Complexity'].isin([10]))]
df_klcluster2['Method'].replace('KlCluster', 'KlCluster 10', inplace=True)

df_klcluster3 = df_filtered[(df_filtered['Method'] == 'KlCluster') & (df_filtered['Complexity'].isin([15]))]
df_klcluster3['Method'].replace('KlCluster', 'KlCluster 15', inplace=True)

df_combined = pd.concat([df_klcluster1, df_klcluster2, df_klcluster3, df_tmm, df_aca, df_haca, df_gmm])


#print(str(df_combined.to_string()))
# Group by TAG and Method, and aggregate by taking the ARITHMETRIC MEAN
df_grouped = df_combined.groupby(['TAG', 'Method'], as_index=False).mean()

# Group by TAG and Method, and aggregate by taking the GEOMETRIC MEAN
#df_grouped = df_combined.groupby(['TAG', 'Method']).agg(lambda x: gmean(x) if x.name not in ['TAG', 'Method'] else x.iloc[0]).reset_index()
#print(df_grouped_)


# df_grouped = df_combined.groupby(['Method'], as_index=False).agg({
#     measure : gmean,
#     "Complexity" : min,
#     "TAG" : max,
# })
#[{measure}, complexity].apply(lambda x: gmean(x, axis=0))")
#df_grouped['accuracy'] = df_grouped.apply(lambda row: gmean(row.dropna()), axis=1)
#np.exp(np.log(df_grouped.prod(axis=1))/df_grouped.notna().sum(1))


print(df_grouped)
# Pivot the table so that each method is a column and the index is TAG
df_pivot = df_grouped.pivot(index='TAG', columns='Method', values=measure)

method_order = ['KlCluster 5', 'KlCluster 10', 'KlCluster 15', 'Tmm', 'Aca', 'Haca', 'Gmm']

# Pivot the table so that each method is a column and the index is 'TAG', and reindex the columns
df_pivot = df_grouped.pivot(index='TAG', columns='Method', values=measure)



# Reorder the columns according to the specified method order
df_pivot = df_pivot[method_order]

# Generate output for each method in the specified format
#print(df_pivot)
print(df_pivot.to_csv())
