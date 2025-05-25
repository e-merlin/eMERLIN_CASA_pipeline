# Configuration

This page explains how to configure the e-MERLIN CASA Pipeline (eMCP) for your specific needs.

## Inputs File

The main configuration for eMCP is through the `inputs.ini` file. Here's a breakdown of the important parameters:

### Basic Settings

```ini
[inputs]
# Project name
inbase          = project_name

# Directory with FITS-IDI files
fits_path       = /path/to/fits/files/

# User email (notifications)
notify          = user@domain.com
```

### Source Lists

Define the sources in your dataset:

```ini
# Main target sources
targets         = M82, M81

# Phase calibrators
phscals         = J0958+6533, J1048+7143

# Flux calibrator
fluxcal         = 1331+305

# Bandpass calibrator
bpcal           = 1407+284

# Point source calibrator (optional)
ptcal           = 

# Reference antenna
refant          = Mk2
```

### Pipeline Flow Control

Control which steps are executed:

```ini
# Enable flagging with AOflagger
flag_aoflagger  = 1

# Automatic a-priori flagging
flag_apriori    = 1

# Manual flagging
flag_manual     = 1

# Time and channel averaging
average         = 1
```

## Configuration Files

### Initialization

The easiest way to get started with configuration files is to use the `init` command:

```bash
# Create configuration files in your current directory
emcp init
```

This will create two important files:
1. `inputs.ini` - For project settings and source lists
2. `default_params.yaml` - For detailed pipeline parameters

### default_params.yaml

More detailed settings are in `default_params.yaml`. This file contains parameters for:

- Flagging thresholds
- Calibration settings
- Imaging parameters
- Weblog options

Key parameters you might want to adjust:

```yaml
flag:
  shadow_limit: 0.0
  quack:
    interval: 5.0
    mode: beg

gaincal_final:
  p_solint: int  # Solution interval for phase calibration
  ap_solint: 60s  # Solution interval for amplitude/phase calibration

first_images:
  wsclean:
    -niter: 10000
    -auto-threshold: 1.0
    -auto-mask: 5.0
```

### YAML Configuration

The YAML format offers several advantages:

- More readable syntax with less visual noise
- Support for comments to document parameters
- Better structure for complex nested configurations

The pipeline looks for configuration files in this order:
1. `./default_params.yaml` (local YAML)
2. Package default YAML (installed with the package)

## Flagging Files

Two main flagging files can be used:

1. `manual.flags` - For pre-averaging flagging
2. `manual_avg.flags` - For post-averaging flagging

Example flag file content:

```bash
mode='manual' antenna='Mk2'
mode='manual' timerange='2019/04/08/12:00:00~12:30:00'
mode='manual' spw='2:0~5'
mode='tfcrop' correlation='ABS_XX,YY' ntime=51.0
```
