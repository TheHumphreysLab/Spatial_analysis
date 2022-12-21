## Python script to compare gene imputation results from SpaGE, Tangram and gimVI
Please installl [SpaGE](https://github.com/tabdelaal/SpaGE), [Tangram](https://github.com/broadinstitute/Tangram) and [scvi-tools](https://github.com/scverse/scvi-tools) before running this script. 
Below are some example codes for using the Leave-One-Out Cross Validation (LOOCV) approach to compare the performance of different spatial gene imputation pipelines. It requires four files as input: 1) ***scRNA cellxgene count***; 2) ***scRNA metadata***; 3) ***spatial cellxgene count***; 4) ***spatial metadata***.  
#### 1. SpaGE:
```bash
# Please change the gene names and working dir in imputation.py
python imputation.py\ 
          --sc_rna /home/data/Qiao/MouseIRI_control.dge.txt.gz\ 
          --sc_rna_meta /home/data/Qiao/GSE139107_MouseIRI.metadata.txt.gz\ 
          --segmentation_cell_stats /home/data/Qiao/segmentation_cell_stats.csv\ 
          --segmentation_counts /home/data/Qiao/segmentation_counts.tsv\
          --output_dir /home/data/Qiao/impute_result\
          --method spaGE
```
#### 2. gimVI:
```bash
# Please change the gene names and working dir in imputation.py
python imputation.py\ 
          --sc_rna /home/data/Qiao/MouseIRI_control.dge.txt.gz\ 
          --sc_rna_meta /home/data/Qiao/GSE139107_MouseIRI.metadata.txt.gz\ 
          --segmentation_cell_stats /home/data/Qiao/segmentation_cell_stats.csv\ 
          --segmentation_counts /home/data/Qiao/segmentation_counts.tsv\
          --output_dir /home/data/Qiao/impute_result\
          --method gimVI
```
#### 3. tangram:
```bash
# Please change the gene names and working dir in imputation.py
python imputation.py\ 
          --sc_rna /home/data/Qiao/MouseIRI_control.dge.txt.gz\ 
          --sc_rna_meta /home/data/Qiao/GSE139107_MouseIRI.metadata.txt.gz\ 
          --segmentation_cell_stats /home/data/Qiao/segmentation_cell_stats.csv\ 
          --segmentation_counts /home/data/Qiao/segmentation_counts.tsv\
          --output_dir /home/data/Qiao/impute_result\
          --method tangram
```
          
