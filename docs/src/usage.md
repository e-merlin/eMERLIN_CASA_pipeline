# Usage

This guide explains how to use the e-MERLIN CASA Pipeline (eMCP) for data processing.

## Quick Start

### Setting Up a New Project

To start a new project, you can use the eMCP pipeline by following these steps:

1. Install the pipeline (see [Installation](installation.md))
2. Create a working directory for your project
3. Initialize a new project with configuration files:

```bash
# Initialize a new project with default configuration files
emcp init

# Force overwrite existing files if needed
emcp init --force
```

This will create:
- `inputs.ini` - The main inputs file where you define project name, data paths, and sources
- `default_params.yaml` - Configuration parameters with default values

### Running the Pipeline

If you have received calibrated data from the observatory and want to refine the calibration:

1. Initialize your project directory with the configuration files
2. Modify `inputs.ini` with your project details and source names
3. [Optionally] Modify `default_params.yaml` or add manual flags to `manual_avg.flags` with your desired values
4. Run the calibration steps:

```bash
# Standard execution with installed package
emcp -r calibration

# Using Docker
docker run -it --rm -v $(pwd):/data emerlin/emcp:latest emcp -r calibration
```

## Command Line Interface

The eMCP package is run with the `emcp` command after installation.

### Basic Usage

```bash
# Show help
emcp -h

# Show version
emcp -v

# List available steps
emcp -l

# Run with custom inputs file
emcp -i my_inputs.ini

# Run specific steps
emcp -r run_importfits,flag_apriori

# Skip specific steps
emcp -s plot_data,flag_target
```

### Running with Parallelization

You can use MPI to process data in parallel:

```bash
# Run with parallelization using MPI
mpicasa -n <num_cores> emcp -r calibration
```

## Available Steps

The pipeline consists of various processing steps that can be selected individually or by categories:

### Pre-processing Steps

- `run_importfits` - Import FITS-IDI files
- `flag_aoflagger` - Run AOflagger for automatic flagging
- `flag_apriori` - Apply a-priori flags
- `flag_manual` - Apply manual flags
- `average` - Time and frequency averaging
- `plot_data` - Generate diagnostic plots
- `save_flags` - Save flag status

### Calibration Steps

- `restore_flags` - Restore flag status
- `flag_manual_avg` - Apply manual flags to averaged data
- `init_models` - Initialize calibrator models
- `bandpass` - Initial bandpass calibration
- `initial_gaincal` - Initial gain calibration
- `fluxscale` - Flux density scaling
- `bandpass_final` - Final bandpass calibration
- `gaincal_final` - Final gain calibration
- `applycal_all` - Apply calibration
- `flag_target` - Flag target sources
- `plot_corrected` - Plot calibrated data
- `first_images` - Make first images
- `split_fields` - Split sources

## Examples

For beginners or typical cases, these command examples will be useful:

1. Run all calibration steps (ideal for observatory-processed data for which you want to tweak the calibration parameters):

```bash
emcp -r calibration
```

2. Run all pipeline steps (you will need the raw FITS-IDI files for the initial step):

```bash
emcp -r all
```

3. Run only pre-processing steps:

```bash
python eMERLIN_CASA_pipeline/eMERLIN_CASA_pipeline.py -r pre_processing
```

## Working with Inputs Files

The pipeline requires an inputs file (default: `inputs.ini`) that specifies:
- Data paths
- Target sources
- Calibrator sources
- Processing options

A template of this file is included in the repository. Copy it to your working directory and modify as needed.
