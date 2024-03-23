# tb-calitas

CALITAS wrapper for Cas9 off-target prediction 

**Install dependencies**\
Clone repository
```
git clone git@github.com:tomebio/tb-calitas.git
cd tb-calitas
```
Create conda environment
```
mamba env create -f environment.yml
```
**WAREHOUSE credentials**\
The pipeline expects WAREHOUSE credentials to be stored in `~/.env`
```
WAREHOUSE_USERNAME="your_username"
WAREHOUSE_PASSWORD="your_passwd"
WAREHOUSE_URL="postgres-warehouse.tome.benchling.com"
```
**Configuration**\
The following are default parameters 
```
params.spacers = "${workflow.projectDir}/test/spacers.txt"
params.pam = 'nrg'
params.fasta = "${workflow.projectDir}/test/chr21.sample.fa"
params.assembly_name = "chr21_sample"
params.species = "human"
params.variants = null
params.max_guide_diffs = 4
params.max_pam_mismatches = 0
params.max_gaps_between_guide_and_pam = 0
```
User-defined parameters can be read from `nextflow.config`
```
params {
       spacers = "./my_spacers.txt"
       fasta = "/data/refs/GRCh38.p14/GRCh38.p14.fasta"
       assembly_name = "GRCh38/hg38"
       species = "human"
}
```
The main input is `my_spacers.txt`, which contains the list of spacer ids (SPXXXX) from benchling, one spacer per line, no header, for example
```
SP2114
SP2115
```
**Run the pipeline**\
Activate conda environment
```
conda activate tb-calitas
```

```
cd /path/to/analysis
nextflow run /path/to/tb-calitas/main.nf [-c nextflow.config]
```
Note: if `nextflow.config` is present in the working directory, it will be loaded automatically