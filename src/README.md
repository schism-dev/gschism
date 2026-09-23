# How to run gschism on TACC VISTA

1. [Command line and arguments](#command-line-and-arguments)
2. [Setting up the Python environment](#setting-up-the-python-environment)
3. [Running it using Slurm's queue](#running-it-using-slurms-queue)

Commands below assume you are in the `src` directory.

## Command line and arguments

```bash
python python/main.py inDir [--outDir OUTDIR] [--nGpus N]
```

| Argument | Required | Meaning |
|---|---|---|
| `inDir` | yes | Path to the SCHISM case directory (`hgrid.gr3`, `hgrid.ll`, `vgrid.in`, `param.nml`, `bctides.in`, …). |
| `--outDir` | no | Output directory. Default: `inDir/gpuoutput`. |
| `--nGpus` | no | Number of GPUs to use (default: `1`). Fails if you ask for more than are available. |


## Setting up the Python environment

gschism needs a set of Python packages before it will run. Most are standard scientific packages.

We use a micromamba environment named `cupy311`. The package list is in `environment.yaml`. Create it once:

```bash
export MAMBA_EXE="$HOME/bin/micromamba"
eval "$("$MAMBA_EXE" shell hook --shell bash)"
micromamba env create -n cupy311 -f environment.yaml
```

If the environment already exists, skip that. The Slurm example below activates it for you. For an interactive shell:

```bash
micromamba activate cupy311
```

## Running it using Slurm's queue

Example batch script: `run_examples/vista_gb.sbatch` (queue `gb`, GH200). It activates `cupy311`, loads the modules it needs, and calls `python/main.py` with the case path.

Before submitting, edit `SRC_DIR` in that file so it is the absolute path to this `src` directory on VISTA (the directory that contains `python/main.py`). The script hardcodes that path because Slurm runs a spool copy of the batch file, so `$0` cannot be used to find the checkout.

```bash
sbatch run_examples/vista_gb.sbatch /scratch/10922/jmglez/tests/RUN02i
sbatch run_examples/vista_gb.sbatch /scratch/10922/jmglez/tests/RUN02i /scratch/10922/jmglez/tests/RUN02i/gpuoutput
```

Logs: `gschism.<jobid>.out` and `gschism.<jobid>.err`. See the sbatch file for modules and library setup.
