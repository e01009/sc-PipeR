--
ver0
config/default.yaml with all user parameters.
Config validation before analysis starts.
Input support:10x directory via Read10X().
Seurat .rds object.

Seurat object creation.
Mitochondrial percentage calculation with organism-aware gene pattern:human: ^MT-
mouse: ^mt-

Cell filtering:min_features
max_features
max_percent_mt

Normalization.
Variable features.
Scaling.
PCA.
UMAP.
Neighbor graph.
Clustering.
Marker detection.
Export:filtered Seurat object .rds
marker table .csv
cell metadata .csv
cluster counts .csv
QC summary before/after filtering .csv
parameters used .yaml
sessioninfo.txt
plain text run summary