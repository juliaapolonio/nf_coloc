# nf_coloc
Nextflow pipeline to extract significant windows and perform colocalization of GWAS and QTL data


## Info and instructions
### Setup
First, download and install [Nextflow](https://docs.seqera.io/nextflow/install) (preferably on Conda) and [Docker](https://www.docker.com/get-started/)/[Singularity](https://docs.sylabs.io/guides/3.0/user-guide/installation.html)

See run_command.sh for a command example to launch the pipeline. Other relevant flags are -resume and -with-tower (check Nextflow docs on [tower](https://docs.seqera.io/platform-enterprise/supported_software/agent/overview) and [resume](https://seqera.io/blog/demystifying-nextflow-resume/)).

### Inputs
The samplesheet is a .tsv (it HAS to be tab-separated) file with filename and path:

```
sample  qtl_file
singlebrain_opc /home/jamorim/data/AD/singlebrain_sceqtl/OPC_for_coloc.tsv
```

GWAS file must have the following columns (including the exact column names):
- snp: rsID of variant
- chr: chromosome (number only)
- bp: base pair location (ensure that gwas and QTL are in the same build, because of window filter)
- freq: minor allele frequency
- b: beta
- se: standard deviation
- p: p-value
- N: number of samples

QTL file must have the same columns and additionally:
- symbol: a gene/protein/metabolite ID (it can be any - symbol, ENSG, etc., but this will be the locuszoom plot title)

The .clumped file is the output of Plink's clump. Do it before running this pipeline (will be included soon for EUR files)
