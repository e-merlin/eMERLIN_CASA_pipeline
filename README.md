# eMERLIN CASA Pipeline

## Contents

1. [Description](#description)
1. [Installation](#installation)
   - [Conda Installation](#conda-installation)
   - [Pip Installation](#pip-installation)
   - [Containers Installation](#containers)
     - [Singularity](#singularity)
     - [Apptainer](#apptainer)
     - [Docker](#docker)
1. [Quick start](#quick-start)
1. [Usage](#usage)
1. [Additional information](#additional-information)
1. [FAQ](#faq)

## Description

The e-MERLIN CASA Pipeline (eMCP) is a python pipeline working on top of [CASA](https://casa.nrao.edu/) to process and calibrate interferometric data from the [e-MERLIN](http://www.e-merlin.ac.uk/) array. Access to data information, statistics and assessment plots on calibration tables and visibilities can be accessed by the pipeline weblog, which is updated in real time as the pipeline job progresses. The output is calibrated data and preliminary lookup images of the relevant fields. It can calibrate mixed mode data that includes narrow-band high spectral resolution spectral windows for spectral lines, and also special observing modes as pseudo-wideband observations. Currently no polarization calibration is performed.

The pipeline uses YAML for all data serialization, including calibration tables, flag statistics, and pipeline information. This provides human-readable data files that can be easily inspected and modified if needed.

## Installation

The e-MERLIN CASA Pipeline (eMCP) requires:
- Python 3.8-3.12. Use Python 3.12 for current modular CASA containers.

Optional requirements:
- aoflagger v2.9+ (needed for L-band data)
- wsclean (alternative to tclean for faster imaging)

CASA does not support all versions of python for all operating systems, so check the documentation in [casadocs](https://casadocs.readthedocs.io/en/stable/notebooks/introduction.html#Compatibility), in particular the Compatibility > Modular CASA section for the right python version for you.

### Method 1: Conda with environment file (Recommended)
Installs Python 3.10 by default and all the dependencies. You still need to install aoflagger and wsclean separately if needed.

We recommend using [mamba](https://mamba.readthedocs.io/en/latest/index.html) for faster installation. If you have conda but not mamba, replace `mamba` by `conda` in the command below:

```bash
git clone https://github.com/e-merlin/eMERLIN_CASA_pipeline.git
cd eMERLIN_CASA_pipeline
git checkout casa6
mamba env create -f environment.yaml -y
conda activate emcp
```

### Method 2: Direct pip installation (no cloning required)

Make sure to have the right python version and pip in your system or environment. You can use virtualenv, conda (`mamba create -n emcp python=3.10 pip -y`) or any other method

You have two alternatives:

2.1 Direct pip installation

```bash
pip install git+https://github.com/e-merlin/eMERLIN_CASA_pipeline.git@casa6
```

2.2 Local installation (allows you to modify the code if needed)

```bash
git clone https://github.com/e-merlin/eMERLIN_CASA_pipeline.git
cd eMERLIN_CASA_pipeline
git checkout casa6
pip install .
```
For development, use `pip install -e .` to install in editable mode.

## Containers

The container includes eMCP, modular CASA, WSClean, and AOFlagger. Run it from
the directory containing `inputs.ini`. CASA data is not baked into the image; the
examples below bind the usual host location, `$HOME/.casa/data`.

### Singularity

```bash
singularity pull emerlin_casa.sif docker://ghcr.io/e-merlin/emerlin_casa_pipeline:base
singularity exec --cleanenv \
  --home "$PWD:/work" \
  --bind "$HOME/.casa/data:/root/.casa/data" \
  emerlin_casa.sif emcp -i inputs.ini -r all
```

### Apptainer

```bash
apptainer pull emerlin_casa.sif docker://ghcr.io/e-merlin/emerlin_casa_pipeline:base
apptainer exec --cleanenv \
  --home "$PWD:/work" \
  --bind "$HOME/.casa/data:/root/.casa/data" \
  emerlin_casa.sif emcp -i inputs.ini -r all
```

### Docker

```bash
docker pull ghcr.io/e-merlin/emerlin_casa_pipeline:base
docker run --rm -it \
  -v "$PWD:/work" \
  -v "$HOME/.casa/data:/root/.casa/data" \
  ghcr.io/e-merlin/emerlin_casa_pipeline:base emcp -i inputs.ini -r all
```

If your CASA data is stored elsewhere, replace `$HOME/.casa/data` with that
path.

If Apptainer/Singularity reports that `squashfuse` is missing and extracts the
SIF to a temporary sandbox, install `squashfuse` on the host system. This cannot
be fixed from inside the container image. As a workaround, set
`APPTAINER_TMPDIR` and `APPTAINER_CACHEDIR` to a filesystem with enough free
space before running `apptainer exec`.

## Quick start

If you have received calibrated data from the observatory and you want to refine the calibration, you can:

1. [Optionally] Modify `default_params.yaml` or add manual flags to `manual_avg.flags` with your desired values.
2. Run one of the following commands depending on your installation method:

- If installed with pip or conda:

  ```bash
  emcp -r calibration
  ```

- If using the repository directly:

  ```bash
  emcp -r calibration
  ```

- If using Docker:

  ```bash
  docker run -it --rm \
    -v "$PWD:/work" \
    -v "$HOME/.casa/data:/root/.casa/data" \
    ghcr.io/e-merlin/emerlin_casa_pipeline:base emcp -r calibration
  ```

## Usage

The eMCP package can be run with the `emcp` command (if installed with pip or conda) or through CASA with the main script.

### Command Line Interface (CLI)

```bash
# Show help
emcp -h

# Show version
emcp -v

# List available pipeline steps
emcp -l

# Multiple steps can be space-separated or comma-separated
emcp -r flag_apriori,flag_manual,average
emcp -r flag_apriori flag_manual average

# Skip specific steps (same comma-separated format)
emcp -s plot_data,save_flags

# Use a custom inputs file
emcp -i my_inputs.ini -r flag_apriori,flag_manual
emcp -v

# List available steps
emcp -l

# Initialize a new project (create default_params.yaml and inputs.ini in current directory)
emcp --init

# Initialize a new project and force overwrite of existing files
emcp --init --force

# Run with custom inputs file
emcp -i my_inputs.ini

# Run specific steps
emcp -r run_importfits,flag_apriori

# Skip specific steps
emcp -s plot_data,flag_target
```

### Starting a new project

```bash
# Create a new directory for your project
mkdir my_project
cd my_project

# Initialize with default configuration files
emcp --init

# Edit inputs.ini and default_params.yaml as needed
```

### Normal pipeline execution

When you have in your working directory the file `inputs.ini`:

```bash
# Standard execution
emcp
```

To run the parallelized version using MPI:

```bash
mpicasa -n <num_cores> emcp
```

For more details, check the [full documentation](docs/index.md).

### Optional arguments

Names in capital need to be set by the user:

```text
  -h, --help                     show this help message and exit


  -i INPUTS_FILE
  --inputs INPUTS_FILE
                                 Inputs file to use. Default is inputs.ini


  -r RUN_STEPS
  --run-steps RUN_STEPS
                                 List of steps to run (space or comma-separated). For example:
                                 "flag_apriori flag_manual average" or "flag_apriori,flag_manual,average"
                                 Also accepts "all", "pre_processing",  "calibration" and "imaging"


  -s SKIP_STEPS
  --skip-steps SKIP_STEPS
                                 List of steps to skip (space or comma-separated)


  -l
  --list-steps                   Show list of available steps and exit
  
  --init                         Initialize a new project directory with config files
  --init --force                 Overwrite existing config files when initializing
```

You can get the list of available steps with:

`emcp -l`

```text
pre_processing
    run_importfits
    flag_aoflagger
    flag_apriori
    flag_manual
    average
    plot_data
    save_flags

calibration
    restore_flags
    flag_manual_avg
    init_models
    bandpass
    initial_gaincal
    fluxscale
    bandpass_final
    gaincal_final
    applycal_all
    flag_target
    plot_corrected
    first_images
    split_fields
```

Selection options are any combination of individual step names, `pre_processing`, `calibration`, or `all`.

### Examples of step selection

Specify which steps of the pipeline to run. For example:

1. Run all calibration steps. This is the usual choice for observatory-processed data when you want to refine calibration parameters:

`emcp -r calibration`

1. Run all pipeline steps (you will need the raw FITS-IDI files for the initial step):

   ```bash
   emcp -r all
   ```

1. Run only the pre-processing steps. These are usually executed by the observatory; otherwise you need the raw FITS-IDI files:

   ```bash
   emcp -r pre_processing
   ```

1. Any combination of the steps above, for example:

   ```bash
   emcp -r plot_corrected first_images split_fields
   ```

1. Run all calibration steps except plot_corrected:

   ```bash
   emcp -r calibration -s plot_corrected
   ```

### Running the pipeline interactively

To execute the pipeline from a running CASA instance you need to write in the CASA shell:

```python
run_in_casa = True
pipeline_path = '/path/to/pipeline_path/'   # You need to define this variable explicitly
execfile(pipeline_path + 'eMERLIN_CASA_pipeline.py')
eMCP = run_pipeline(run_steps=['calibration'])
```

Function `run_pipeline` parameters and defaults are: `run_pipeline(inputs_file='./inputs.ini', run_steps=[], skip_steps=[])`. Variables run_steps and skip_steps are python lists of steps as explained above.

## Additional information

### Data Storage

The pipeline uses YAML files for storing all serialized data, including:

- Pipeline information (`eMCP_info.yaml`)
- Calibration tables (`caltables.yaml`)
- Flag statistics (`flagstats_*.yaml`)
- Flux calibration results (`calfluxes.yaml`)

These YAML files are human-readable and can be inspected or modified with any text editor. This format allows for easier debugging and customization compared to binary formats.

- [Documentation](#additional-information)
- [Wiki pages](https://github.com/e-merlin/eMERLIN_CASA_pipeline/wiki)

## FAQ

### How do I open the weblog?

The weblog consist of a series of html files. From the working directory you can open the file `./weblog/index.html` with your preferred web browser.

**How do I know what has been executed?**

You can visit the tab `Pipeline info` in the weblog, where you will find which steps were executed. You will also find a link to the Pipeline log, the CASA log and two files with all the parameters used during the data processing.

**I want to re-run the pipeline to improve the calibration, what do I change?**

There are two main blocks: pre-processing and calibration. Most probably you will only need to repeat the calibration part. Recommended course of action:

- Identify changes you want to include in the data reduction, like changing calibration parameters or adding manual flags.
- Add or edit file `manual_avg.flags` with your flag commands (follow the CASA syntax).
- Edit the file `inputs.ini` if you need to change the sources used or they intend.
- Edit the file `default_params.yaml` changing any parameter the pipeline is using, if needed.
- Run the calibration block of the pipeline with the command:

`emcp -r calibration`

**Which flag files does the pipeline accept and what is the right syntax?**

There are four different flag files accepted by the pipeline:

| Flag file     | Used by step  | Notes |
| ------------- |:-------------:| -----|
| observatory.flags | flag_apriori | Created by the observatory with antenna slewing or other major faults. Please donot edit it yourself. |
| manual.flags | flag_manual | This is meant to flag the unaveraged data set during the pre-processing stage |
| manual_avg.flags | flag_manual_avg | This is meant to flag the averaged data set during the calibration stage |
| manual_narrow.flags | flag_manual_avg | Use this to add flag commands for narrow-band spectral line data set|

For the syntax needed for CASA follow [Basic Syntax Rules](https://casa.nrao.edu/casadocs/casa-5.5.0/global-task-list/task_flagdata/about) in the CASA docs flagdata (end of the section). The main rules are:

1. Use only ONE white space to separate the parameters (no commas). Each key should only appear once on a given command line/string.
2. There is an implicit mode for each command, with the default being 'manual' if not given.
3. Comment lines can start with '#' and will be ignored. The parser used in flagdata will check each parameter name and type and exit with an error if the parameter is not a valid flagdata parameter or of a wrong type.

Example for e-MERLIN:

```text
mode='manual' field='1331+305' antenna='' timerange='10:00:00~10:11:30'
mode='manual' field='' antenna='' timerange='' spw='0:0~30'
mode='manual' field='' antenna='Mk2' timerange='09:05:00~16:27:00'
mode='manual' field='1258-2219' antenna='' timerange='12:57:01~12:59:59'
mode='quack' field='1258-2219,1309-2322' quackinterval=24.
```

Example from the CASA docs:

```python
scan='1~3' mode='manual'
# this line will be ignored
spw='9' mode='tfcrop' correlation='ABS_XX,YY' ntime=51.0
mode='extend' extendpols=True
scan='1~3,10~12' mode='quack' quackinterval=1.0
```

**How do I fill the source names in inputs.ini if I don't know which fields were observed?**

By default you should have all the information from the observatory. But if you only have the FITS-IDI and don't know the source names, you can run the first pipeline step alone `emcp -r run_importfits`. When the execution is finished, open the weblog and go to the tab `Observation Summary` where you will find the fields included in the MS and the listobs file with all the scans.

As a general rule, an observation will have 1331+3030 (3C286) as flux scale calibrator, 1407+2827 (OQ208) as bandpass calibrator and 0319+4130 (3C84) as bright ptcal calibrator. To distinguish between target and phasecal, you should look for alternating scans, and the target is usually the one with longer scans.
