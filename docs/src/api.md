# API Documentation

This page provides documentation for the key modules and functions in the e-MERLIN CASA Pipeline that may be useful for users who want to build custom workflows.

## Main Modules

The eMCP package is organized into several modules:

### eMCP_functions

Contains the core pipeline functions for data processing and calibration.

```python
from eMCP.functions import eMCP_functions as em
```

### eMCP_utils

Contains utility functions for file handling, configuration, and environment setup.

```python
from eMCP.utils import eMCP_utils as emutils
```

### eMCP_plots

Contains functions for generating plots and visualizations.

```python
from eMCP.plots import eMCP_plots as emplt
```

## Key Functions

### Data Import and Setup

```python
# Import FITS-IDI files
em.import_eMERLIN_fitsIDI(eMCP)

# Get metadata from MS file
em.get_msinfo(eMCP, msfile)

# Create directory structure
emutils.create_dir_structure()
```

### Flagging Functions

```python
# Run AOflagger
em.run_aoflagger_fields(eMCP)

# Apply a-priori flags
em.flagdata1_apriori(eMCP)

# Apply manual flags
em.flagdata_manual(eMCP, run_name='flag_manual')

# Save flag status
em.saveflagstatus(eMCP)

# Restore flag status
em.restoreflagstatus(eMCP)
```

### Calibration Functions

```python
# Initialize calibration tables dictionary
em.initialize_cal_dict(eMCP)

# Initialize calibrator models
em.run_initialize_models(eMCP)

# Apply bandpass calibration
em.initial_bp_cal(eMCP, caltables)

# Apply gain calibration
em.initial_gaincal(eMCP, caltables)

# Apply flux scale
em.eM_fluxscale(eMCP, caltables)

# Apply final bandpass calibration
em.bandpass_final(eMCP, caltables)

# Apply final gain calibration
em.gaincal_final(eMCP, caltables)

# Apply all calibration
em.applycal_all(eMCP, caltables)
```

### Plotting Functions

```python
# Plot elevation and UV coverage
em.plot_elev_uvcov(eMCP)

# Generate standard diagnostic plots
emplt.make_4plots(eMCP, datacolumn='data')
```

### Imaging Functions

```python
# Create first images
em.run_first_images(eMCP)
```

## Pipeline Configuration

The pipeline is configured using a dictionary structure that contains all settings and state information:

```python
# Example of creating the eMCP dictionary
eMCP = emutils.start_eMCP_dict(info_dir)

# Load default parameters
defaults_file = './default_params.json'
eMCP['defaults'] = json.loads(open(defaults_file).read())

# Load inputs
inputs = emutils.read_inputs(inputs_file)
eMCP['inputs'] = inputs
```

## Example Custom Script

Here's an example of a simple custom script that uses the eMCP functions:

```python
#!/usr/bin/env python
"""
Custom script to process e-MERLIN data
"""
import os
import json
from eMCP.functions import eMCP_functions as em
from eMCP.utils import eMCP_utils as emutils
from eMCP.plots import eMCP_plots as emplt

# Create directory structure
calib_dir, info_dir = emutils.create_dir_structure()

# Initialize eMCP dictionary
eMCP = emutils.start_eMCP_dict(info_dir)

# Load default parameters
eMCP['defaults'] = json.loads(open('./default_params.json').read())

# Load inputs
eMCP['inputs'] = emutils.read_inputs('./inputs.ini')

# Get MS file info
msfile = eMCP['inputs']['inbase'] + '.ms'
eMCP, msinfo, msfile = em.get_msinfo(eMCP, msfile)

# Run flagging
eMCP = em.run_aoflagger_fields(eMCP)
eMCP = em.flagdata1_apriori(eMCP)

# Plot raw data
em.plot_elev_uvcov(eMCP)
eMCP = emplt.make_4plots(eMCP, datacolumn='data')

print("Custom processing completed!")
```
