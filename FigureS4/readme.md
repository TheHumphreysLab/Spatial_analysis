## Python script to compare gene imputation results from SpaGE, Tangram and gimVI
Please installl [SpaGE](https://github.com/tabdelaal/SpaGE), [Tangram](https://github.com/broadinstitute/Tangram) and [scvi-tools](https://github.com/scverse/scvi-tools) before running this script. 
Below are some example codes for using the Leave-One-Out Cross Validation (LOOCV) approach to compare the performance of different spatial gene imputation pipelines. It requires four files as input:1) ***scRNA cell***  $\times$ ***gene count***; 2) ***scRNA metadata***; 3) ***spatial cell*** $\times$ ***gene count***; 4) ***spatial metadata***.  
#### 1. SpaGE:
```bash
# Please change the gene names and working dir in imputation.py
python imputation.py\ 
          --sc_rna /home/data/MouseIRI_control.dge.txt.gz\ 
          --sc_rna_meta /home/data/GSE139107_MouseIRI.metadata.txt.gz\ 
          --segmentation_cell_stats /home/data/segmentation_cell_stats.csv\ 
          --segmentation_counts /home/data/segmentation_counts.tsv\
          --output_dir /home/data/impute_result\
          --method spaGE
```
#### 2. gimVI:
```bash
# Please change the gene names and working dir in imputation.py
python imputation.py\ 
          --sc_rna /home/data/MouseIRI_control.dge.txt.gz\ 
          --sc_rna_meta /home/data/GSE139107_MouseIRI.metadata.txt.gz\ 
          --segmentation_cell_stats /home/data/segmentation_cell_stats.csv\ 
          --segmentation_counts /home/data/segmentation_counts.tsv\
          --output_dir /home/data/impute_result\
          --method gimVI
```
#### 3. tangram:
```bash
# Please change the gene names and working dir in imputation.py
python imputation.py\ 
          --sc_rna /home/data/MouseIRI_control.dge.txt.gz\ 
          --sc_rna_meta /home/data/GSE139107_MouseIRI.metadata.txt.gz\ 
          --segmentation_cell_stats /home/data/segmentation_cell_stats.csv\ 
          --segmentation_counts /home/data/segmentation_counts.tsv\
          --output_dir /home/data/impute_result\
          --method tangram
```
          
