# Pipeline Steps

This page documents the individual processing steps in the e-MERLIN CASA Pipeline.

## Overview

The pipeline consists of two main phases:

1. **Pre-processing** - Import and prepare data for calibration
2. **Calibration** - Apply calibration routines and generate output products

## Pre-processing Steps

### run_importfits

Imports FITS-IDI files into a CASA Measurement Set (MS).

**Function**: Concatenates all FITS-IDI files in the specified folder into a single MS.

**Inputs**: FITS-IDI files path specified in inputs.ini

**Outputs**: `<inbase>.ms` file and `<inbase>.ms.listobs.txt`

**CASA tasks used**: importfitsidi, fixvis, flagdata (autocorr)

### flag_aoflagger

Runs the AOflagger tool for automated RFI detection and flagging.

**Function**: Uses the external AOflagger tool with optimized strategies for e-MERLIN data.

**Requires**: AOflagger v2.9+ installed on the system

**Notes**: Especially important for L-band data where RFI can be significant

### flag_apriori

Applies a priori flagging based on known issues.

**Function**: Flags autocorrelations, shadowed antennas, edge channels, Lovell-Mk2 baseline issues, etc.

**Default flags**:

- Autocorrelations
- Zero amplitude data
- Shadowed antennas with limit = 0.0
- 5 channels at band edges for each spectral window
- Specific channels known to have issues

### flag_manual

Applies manual flags specified in the `manual.flags` file.

**Function**: Reads the flag commands from the file and applies them to the MS.

**Example format**:

```bash
mode='manual' antenna='Mk2'
mode='manual' timerange='2019/04/08/12:00:00~12:30:00'
mode='manual' spw='2:0~5'
```

### average

Performs time and frequency averaging to reduce data volume.

**Function**: Uses CASA's mstransform task to average data in time and frequency domains.

**Outputs**: `<inbase>_avg.ms` or `<inbase>_avg.mms` if parallelization is enabled

**Parameters controlled in default_params.json**:

- Time averaging interval
- Channel averaging width

### plot_data

Generates diagnostic plots for the raw data.

**Function**: Creates amplitude/phase vs time/frequency plots for all sources.

**Outputs**: PNG files in the weblog directory

### save_flags

Saves the current flag state to allow restoration later.

**Function**: Uses CASA's flagmanager to save flags with versionname 'before_cal_flags'

## Calibration Steps

### restore_flags

Restores flags to the state saved with save_flags step.

**Function**: Uses flagmanager to recover previous flag state.

### flag_manual_avg

Applies manual flags to the averaged data from the `manual_avg.flags` file.

**Function**: Similar to flag_manual but operates on the averaged dataset.

### init_models

Initializes calibrator models for flux and bandpass calibrators.

**Function**: Sets appropriate models for standard calibrators (e.g., 3C286, 3C48).

**Notes**: Uses built-in CASA models or sets point source models with known flux densities.

### bandpass

Performs initial bandpass calibration.

**Function**: Creates initial bandpass calibration table for spectral response correction.

**Outputs**: Initial bandpass calibration table

**CASA tasks used**: bandpass

### initial_gaincal

Creates initial gain calibration tables.

**Function**: Derives delay, phase, and amplitude calibration solutions.

**Outputs**:

- Delay calibration table (K)
- Phase calibration table (G)
- Amplitude calibration table (G)

### fluxscale

Bootstraps flux scale from primary calibrator to secondary calibrators.

**Function**: Uses fluxscale task to transfer flux density scale.

**Outputs**: Flux-scaled gain table and flux density values for secondary calibrators

### bandpass_final

Creates final bandpass calibration with spectral index information.

**Function**: Improves bandpass solution using spectral index information from fluxscale.

**Outputs**: Final bandpass calibration table

### gaincal_final

Final gain calibration including spectral information.

**Function**: Creates final phase and amplitude calibration solutions.

**Outputs**: Final gain calibration tables

### applycal_all

Applies all calibration tables to the dataset.

**Function**: Uses applycal to apply calibration to all sources.

**Outputs**: Fully calibrated data in the CORRECTED_DATA column

### flag_target

Applies automated flagging to calibrated target data.

**Function**: Uses rflag or tfcrop algorithm on CORRECTED_DATA.

### plot_corrected

Generates diagnostic plots for the calibrated data.

**Function**: Creates amplitude/phase vs time/frequency plots for calibrated data.

**Outputs**: PNG files in the weblog directory

### first_images

Creates preliminary images of calibrators and targets.

**Function**: Uses CASA's tclean to make simple images.

**Outputs**: FITS images and PNG previews in the images directory

### split_fields

Splits individual sources into separate Measurement Sets.

**Function**: Uses split task to create individual files per source.

**Outputs**: `<source_name>.ms` files for each selected source
