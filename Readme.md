# tb-calitas

CALITAS wrapper for Cas9 off-target prediction 

**Install dependencies**
- Clone repository
```
git clone git@github.com:tomebio/tb-calitas.git
cd tb-calitas
```
- [Install conda][conda-link]
- Install mamba

```
conda install -c conda-forge mamba
```
- Create conda environment
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
**Input**
- `spacers`: .txt file containing the list of spacer ids (SPXXXX) from benchling, one spacer per line, no header
```
SP2114
SP2115
```
- `fasta`: reference FASTA file, requires both index and dictionary, which can be generated using samtools
```
samtools faidx ref.fa
samtools dict -a <assembly-name> -s <species> -o ref.fa.dict ref.fa
```
The following are default parameters 
```
params.spacers = "${workflow.projectDir}/test/spacers.txt"
params.fasta = "${workflow.projectDir}/test/chr21.sample.fa"
params.pam = 'nrg'
params.max_guide_diffs = 4
params.max_pam_mismatches = 0
params.max_gaps_between_guide_and_pam = 0
params.variants = null
```
User-defined parameters can be read from `nextflow.config`.\
```
params {
       spacers = "my_spacers.txt"
       fasta = "/data/refs/GRCh38.p14/GRCh38.p14.fasta"
}
```
Note:`fasta` must be a full path, otherwise nextflow fails to find a relative path.

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

[conda-link]: https://docs.conda.io/projects/conda/en/latest/user-guide/install/