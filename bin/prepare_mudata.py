import numpy as np
import pandas as pd
import anndata as ad
import mudata as mu
import muon
import sys


mudata_object_path = sys.argv[1]
output_mudata_path = sys.argv[2]



# Read in mudata object
mdata = muon.read(mudata_object_path)


# Select omics data to use for SNF
#mdata = mdata['acylcarnitines', 'fattyacids']


# Save mudata object
mdata.write(output_mudata_path)
