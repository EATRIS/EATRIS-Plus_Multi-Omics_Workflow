import sys
import numpy as np
import pandas as pd
import matplotlib.pylab as plt
from sklearn.cluster import spectral_clustering
from sklearn.metrics import v_measure_score
import snf
from snf import metrics
import seaborn as sns
import scipy.stats as stats
import statsmodels.api as sm
import statsmodels
import muon
import itertools
import os


# Arguments
mudata_object = sys.argv[1]
omics_names_file = sys.argv[2]
silhouette_cutoff = sys.argv[3]


# Read MuData object
mdata = muon.read(mudata_object)

# Read in omics combinations names
omics_names_str = str(omics_names_file)[2:-6]
omics_names = omics_names_str.split("', '")



# Functions used for SNF

## Cross validation
#grid_zaff, grid_labels = cv.snf_gridsearch(data_combined)
#mu_cv, k_cv = cv.get_optimal_params(grid_zaff, grid_lables)

def run_snf(data_list):
    # Calculate affinity networks for single omics
    affinity_networks = snf.make_affinity(data_list, metric = 'euclidean', K=6, mu=0.5)

    # Remove infinite values, set to 0 
    for i in affinity_networks:
        i[~np.isfinite(i)] = 0

    # SNF
    fused_network = snf.snf(affinity_networks)
    
    return(fused_network)

def sort_fused_network(fused_network, labels_array):
    # Make Pandas Dataframes
    df = pd.DataFrame(fused_network)
    df_labels = pd.DataFrame(labels_array)
    df_labels.columns = ["Label"]
    # sort label df
    df_labels = df_labels.sort_values(by=['Label'])
    # sort fused network df with sorted labels
    df = df.reindex(df_labels.index)
    df = df[df_labels.index]
    array = df.to_numpy()
    np.fill_diagonal(array, 0)
    return(array)



# Find common samples in omics combination
obs_list = []
for i in omics_names:
    obs_list.append(list(mdata[i].obs.index))

common_obs = list(set.intersection(*map(set, obs_list)))

# Subset MuData object on common IDs
muon.pp.filter_obs(mdata, common_obs)




# Make list of arrays of the omics combination
data_arrays = []
for i in omics_names:
    data_arrays.append(mdata[i].X)
    
# Run SNF
fused_network = run_snf(data_arrays)

# Find best and second optimal number of clusters
best, second = snf.get_n_clusters(fused_network)

# Perform spectral clustering on the fused network
labels = spectral_clustering(fused_network, n_clusters=best)
labels_second = spectral_clustering(fused_network, n_clusters=second)
    
# Silhouette score
np.fill_diagonal(fused_network, 0)
sil = metrics.silhouette_score(fused_network, labels)
sil2 = metrics.silhouette_score(fused_network, labels_second)

         

# Write silhouette score and name of omics types
with open( str(omics_names) + '_sil_scores.tsv', 'w') as f:
    f.write(str(omics_names))
    f.write('\t')
    f.write(str(len(omics_names)))
    f.write('\t')
    f.write(str(sil))
    f.write('\t')
    f.write(str(sil2))
    f.write('\n')
    
if sil > float(silhouette_cutoff):
    fused_df = pd.DataFrame(fused_network)
    fused_df.to_csv(str(omics_names) + '_fused_network.csv', header=False, index=False)