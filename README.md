
# 'build annotations' pipeline

## description

Given a FASTA file, transcript annotation, and ORF annotation, this pipeline creates several annotation files (genic, convergent, divergent, and intergenic regions) which are useful for other pipelines. It also creates a coverage file of GC% in an 11bp running window. This pipeline is usually called as a subworkflow of other pipelines rather than being run independently. This allows multiple pipelines (e.g. for different assays) to share the same 'build-annotations' pipeline.

## requirements:

### software

- Unix-like operating system (tested on CentOS 7.2.1511)
- Git
- [conda](https://conda.io/docs/user-guide/install/index.html)

### files

- FASTA file of genome

- [BED6](https://genome.ucsc.edu/FAQ/FAQformat.html#format1) format annotation files:
    - ORF annotation
    - transcript annotation

## instructions

**0**. Clone this repository

```bash
git clone https://github.com/winston-lab/build-annotations.git
```

**1**. Make a copy of the configuration file template `config_template.yaml` called `config.yaml`, and edit `config.yaml` to suit your needs.

```bash
# navigate into the pipeline directory
cd build-annotations

# create the configuation file
cp config_template.yaml config.yaml

vim config.yaml    # or use your favorite editor

# make your edits
```

**2**. **STOP!** If you only need this pipeline as a subworkflow of another pipeline, you're done! If you really want to run this pipeline independently, create and activate the `build_annotations` virtual environment for the pipeline using conda. The virtual environment creation can take a while.

```bash
# create the snakemake_default environment
conda env create -v -f envs/build-annotations.yaml

# activate the environment
source activate build_annotations

# to deactivate the environment
# source deactivate
```

**3**. Do a dry run of the pipeline to see what files will be created.

```bash
snakemake -p --dry-run
```

**4**. Run the pipeline using the above command, omitting the `--dry-run` flag. 

```bash
snakemake -p
```

