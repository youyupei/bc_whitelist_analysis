This folder contains some script/notebook for the analysis in BLAZE paper

### Code used for running Cell Ranger，BLAZE and Sockeye to get the whitelist:
#### Cell Ranger
```
[to be added]
```
#### BLAZE
```
python3 get_raw_bc.py --expect-cells=500 --threads=32 <path/to/fastq_pass>
```
#### Sockeye
```
snakemake --use-conda --configfile config/config.yml -pr all
```
**Note**: We used all the default parameter in the config file from Sockeye.

### Running FLAMES
Code used for running FLAMES with specific whitelist is in `FLAMES.r`
