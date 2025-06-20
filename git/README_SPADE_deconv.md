
### SpaGCN.py 

[SpaGCN github repository](https://github.com/jianhuupenn/SpaGCN)
##### Input : 
- preprocessed Spatial expression matrix, .csv
- spots locations, .csv
- histology img, .tiff
##### Output : 
- Spatial expression matrix with additional column = cluster annotation for each spot, results.h5ad
- Cluster visualisation, pred.png

-----
### deconv_SPADE.Rmd

[SPADE github repository](https://github.com/YyLu5/SPADE)
##### Input : 
- preprocessed + clustered spatial expression matrix, .rds
- spots locations, .csv
- single cell reference matrix, .rds
- Lists of marker genes (one list per cell type), .rds
##### output : 
- Deconvolution matrix (spots X cell types), .csv
- Pie plot visualisation, .pdf
- Co-localisation matrix, .pdf

