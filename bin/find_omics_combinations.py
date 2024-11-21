import numpy as np
import sys
import os
import pandas as pd
import muon
import itertools


# Parameters
mudata_object_path = sys.argv[1]
min_n_combi = int(sys.argv[2])
max_n_combi = int(sys.argv[3])

# Read in mudata
mdata = muon.read(str(mudata_object_path))

# Get all combinations of omics data into lists of numpy arrays
possible_combinations = []
for i in range(min_n_combi, max_n_combi + 1):
    possible_combinations.append(list(itertools.combinations(mdata.mod.keys(), i)))

# Save combinations of omics data to simple txt files
for i in possible_combinations:
    for j in i:
        filename = str(j) + ".txt"
        with open (filename, "w") as f:
            f.write(str(j))