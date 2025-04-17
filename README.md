This repo is designed to compare a simple workflow
between snakemake and nextflow.

First, download and move into the GitHub repository.

```
git clone https://github.com/ccpowers-NOAA/PEMAD-PBB-fastp.git
cd PEMAD-PBB-fastp
```

The raw data are downloaded as part of the repository (`00_reads`).

## Snakemake 

To run snakemake, first setup your snakemake
conda environment

`mamba env create -f envs/snakemake-8.20.3.yaml -n snakemake-8.20.3`

and install the following dependency

`pip install snakemake-executor-plugin-slurm --user`

Execute the pipeline interactively by running
```
# Before running anything interactively, get off the head node!
srun -p standard -c 2 --mem=5GB --pty /bin/bash
mamba activate snakemake-8.20.3.yaml
snakemake -j1 --use-conda
```

or in the queue by running

`sbatch snake_manager.sh`

## Nextflow

To run nextflow, we can use SEDNA's nextflow
conda environment

```
# Before running anything interactively, get off the head node!
srun -p standard -c 2 --mem=5GB --pty /bin/bash
mamba activate nextflow-24.04.4
```

Execute the pipeline interactively by running

`nextflow -c nf.config.NOAA_SEDNA run main.nf -entry QC -with-dag flowchart.png`

or in the queue by running

`sbatch nextflow_job.sh`

We have provided two Nextflow config files. `nf.config.NOAA_SEDNA` will submit jobs
to the SLURM scheduler, while `nf.config.local` will only use local compute resources
(i.e. the CPUs and memory in an interactive session).

In the config file, there are two ways to make Nextflow use conda. You can either build a new, temporary conda environment from our provided `fastp-0.24.0.yaml`. This may take several minutes.

```
conda = "${projectDir}/envs/fastp-0.24.0.yaml"
```

Or you can have Nextflow use an existing conda environment.

```
conda = '/home/dmacguigan/.conda/envs/fastp-0.24.0'
```