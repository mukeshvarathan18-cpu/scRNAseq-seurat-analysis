# convert_h5ad_to_10x.py
# Convert .h5ad single-cell dataset to 10X format

import scanpy as sc

# load dataset
adata = sc.read_h5ad("dataset.h5ad")

# write 10X files
sc.write_10x_mtx(
    "converted_10x",
    adata,
    var_names="gene_symbols",
    overwrite=True
)

print("Conversion completed")