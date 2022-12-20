This folder contains some script/notebook for the analysis in BLAZE paper. The jupyter notebooks (`.ipynb`) contain the comparison among the barcode whitelists from Cell Ranger (using short reads), BLAZE and Sockeye. All the whitelists were used as input to FLAME and the R Markdown(`.Rmd`) files contain the code for generating the UMAP plots.

### Code used for running Cell Ranger，BLAZE and Sockeye to get the whitelist:
#### Cell Ranger
```
see https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/tutorial_ct 
```
#### BLAZE
```
python3 blaze.py --expect-cells=500 --threads=32 <path/to/fastq_pass>
```
#### Sockeye
```
snakemake --use-conda --configfile config/config.yml -pr all
```
**Note**: We used all the default parameter in the config file from Sockeye.

### Running FLAMES
Code used for running FLAMES with specific whitelist is in `FLAMES.r`
