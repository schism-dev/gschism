# How to run gschism on TACC VISTA

This file explains how to set up gschism on VISTA once, and then how to run a SCHISM case on VISTA's GPU nodes.

1. [Getting the code](#getting-the-code)
2. [Setting up the Python environment](#setting-up-the-python-environment)
3. [Running a case with Slurm](#running-a-case-with-slurm)
4. [Command line and arguments](#command-line-and-arguments)
5. [Appendix: Python packages](#appendix-python-packages)



## Getting the code

On a VISTA login node, clone the repository into your home directory:

```bash
mkdir -p "$HOME/src"
cd "$HOME/src"
git clone https://github.com/schism-dev/gschism.git
cd gschism/src
```

All commands below assume you are in this `src` directory, `$HOME/src/gschism/src`. The first run writes compiled files here, so keep the checkout somewhere you can write to.



## Setting up the Python environment

gschism works with Python 3.11 plus a fixed set of packages, each at a specific version. The list is in `environment.yaml`. The [appendix](#appendix-python-packages) describes each package.

The usual way to deal with packages and Python versions is an environment. A Python environment is a directory with its own Python and its own packages. While it is active, `python` means that Python, and the system Python and other environments are left alone.

Among the environment managers (Conda, Mamba, Micromamba), we use Micromamba because it is a single executable file: copying it into your home directory is the whole installation. Conda and Mamba are full installations with a base Python of their own. All three read the same `environment.yaml`.

The steps below are done once, on a VISTA login node.

### How to install Micromamba

VISTA has ARM processors, so download the ARM (`linux-aarch64`) build. These lines put the program at `$HOME/bin/micromamba`, which is where the batch script expects it:

```bash
mkdir -p "$HOME/bin"
curl -Ls https://micro.mamba.pm/api/micromamba/linux-aarch64/latest | tar -xj -C "$HOME" bin/micromamba
```

Then set up your shell for it:

```bash
"$HOME/bin/micromamba" shell init --shell bash --root-prefix "$HOME/.local/share/mamba"
source ~/.bashrc
```

This adds a few lines to `~/.bashrc`. Every new shell then has the `micromamba` command, and environments are kept under `$HOME/.local/share/mamba`, which is where the batch script looks for them. A Slurm job starts with the settings of the shell it was submitted from, so the batch script finds them as well.

### Create the environment

```bash
micromamba env create -n cupy311 -f environment.yaml
```

This downloads a few GB into `$HOME`, which has a disk quota. `quota` shows how much space you have left.

If the install stops while building `pnetcdf`, delete that line from `environment.yaml` and run the command again. gschism does not use `pnetcdf`.

The batch script activates the environment by itself. To use it in a shell of your own:

```bash
micromamba activate cupy311
```



## Running a case with Slurm

On VISTA, programs that need a GPU run on compute nodes, through Slurm, the queue system. `run_examples/vista_gb.sbatch` asks for one GH200 node in the `gb` queue for 20 minutes. It loads the modules gschism needs (GCC, CUDA 13.1, METIS), activates `cupy311`, and runs `python/main.py` on the case you give it.

Before the first run, open the script and set `SRC_DIR` to the full path of this `src` directory. The script cannot work it out by itself, because Slurm runs a copy of the script from a different directory. If the case needs more than 20 minutes, change the `-t` line.

Submit a case like this. The second argument is the output directory, and can be left out:

```bash
sbatch run_examples/vista_gb.sbatch $SCRATCH/tests/RUN02i
sbatch run_examples/vista_gb.sbatch $SCRATCH/tests/RUN02i $SCRATCH/tests/RUN02i/gpuoutput
```

The job writes its screen output to `gschism.<jobid>.out` and its errors to `gschism.<jobid>.err`, in the directory you submitted from.

The first run also compiles the C++ parts of gschism and keeps the results inside `src`. Later runs reuse them.



## Command line and arguments

This is the command the batch script runs:

```bash
python python/main.py inDir [--outDir OUTDIR] [--nGpus N]
```

| Argument   | Required | Meaning                                                                                                |
| ---------- | -------- | ------------------------------------------------------------------------------------------------------ |
| `inDir`    | yes      | Path to the SCHISM case directory (`hgrid.gr3`, `hgrid.ll`, `vgrid.in`, `param.nml`, `bctides.in`, …). |
| `--outDir` | no       | Output directory. Default: `inDir/gpuoutput`.                                                          |
| `--nGpus`  | no       | Number of GPUs to use (default: `1`). Fails if you ask for more than are available.                    |

It needs a GPU and the modules that the batch script loads, so on VISTA use the batch script instead of running it on a login node.



## Appendix: Python packages

What each package in `environment.yaml` is, how common it is, and where gschism uses it.

| Package | What it is | How common | Where gschism uses it |
| --- | --- | --- | --- |
| `python` 3.11 | The Python interpreter. | A common version for scientific work. | Everywhere. |
| `numpy` | Arrays and numerical routines on the CPU. | Standard. | Throughout, for arrays on the CPU. |
| `scipy` | Scientific routines built on NumPy. | Standard. | `scipy.spatial.Delaunay`, when reading atmospheric forcing (sflux). |
| `matplotlib` | Plotting. | Standard. | Loaded at startup by `python/main.py` and `debug/monitor.py`, so it is needed even when nothing is plotted. |
| `netcdf4` | Reads and writes NetCDF files. | Standard in ocean and atmosphere modeling. | Reading sflux and writing model output. |
| `cftime` | Dates in the calendars used by climate and ocean data. | Standard alongside `netcdf4`. | Not directly. `netcdf4` needs it. |
| `cupy-cuda13x` | NumPy-like arrays on an NVIDIA GPU, built for CUDA 13. | The most common GPU array library for Python. | The model arrays, the kernels, and the solver. |
| `fastrlock` | A fast lock. | Rarely installed on its own. | Not directly. CuPy needs it. |
| `pybind11` | Connects C++ code to Python. | The usual tool for this. | Compiling the grid reader (`grid/read_and_process_grid.cpp`) and the solver loop. |
| `f90nml` | Reads Fortran namelist files. | Common around Fortran models such as SCHISM, rare elsewhere. | Reading `param.nml`. |
| `line_profiler` | Measures time spent on each line of code. | A development tool. | Loaded at startup by `python/main.py`, so gschism will not start without it. |
| `pip` | Python's own package installer. | Standard. | Only while the environment is created, to install the packages listed under `pip:` in `environment.yaml`. |
| `Cython` | Compiles Python-like code to C. | A standard build tool. | Not used. |
| `mpi4py` | MPI (message passing between processes) for Python. | Standard on clusters. | Not used. |
| `pnetcdf` | Python interface to the PnetCDF parallel I/O library. | Uncommon. It is built from source during installation and needs the PnetCDF C library. | Not used. |
| `typing_extensions` | Newer type-hint features for older Python versions. | A standard utility. | Not used. |
