# Troubleshooting

This page provides solutions to common issues encountered when using eMCP.

## Common Issues

### Installation Problems

#### Missing Dependencies

**Problem**: Error messages about missing packages when running the pipeline.

**Solution**: Make sure you've installed all dependencies. Using conda is recommended:

```bash
conda env create -f environment.yml
conda activate emcp
```

#### CASA Module Import Errors

**Problem**: Errors like "No module named casatasks".

**Solution**: Make sure you're running in the correct environment:

```bash
conda activate emcp
python -c "import casatasks"  # Test if CASA is correctly installed
```

### Runtime Issues

#### AOflagger Not Found

**Problem**: Error message "AOflagger executable not found in PATH".

**Solution**: 

1. Install AOflagger:

   ```bash
   conda install -c conda-forge aoflagger
   ```

2. Or edit default_params.json to disable it:

   ```json
   {
     "run": {
       "aoflagger": false
     }
   }
   ```

#### Memory Issues

**Problem**: Pipeline crashes with memory errors.

**Solution**:

1. Process smaller chunks of data at a time:

   ```bash
   # Run pre-processing steps only
   emcp -r pre_processing
   
   # Then run calibration steps
   emcp -r calibration
   ```

2. Use more time/channel averaging in default_params.json:

   ```json
   {
     "averaging": {
       "timeavg": 8,
       "chanavg": 8
     }
   }
   ```

### Calibration Issues

#### Bad Bandpass Solutions

**Problem**: Bandpass calibration shows extreme values or fails.

**Solution**:

1. Check if the bandpass calibrator has enough SNR
2. Try flagging more aggressively:

   ```bash
   # Create a manual_avg.flags file with tfcrop settings
   echo "mode='tfcrop' correlation='ABS_XX,YY' ntime=51.0" > manual_avg.flags
   
   # Run with manual flagging before bandpass
   emcp -r restore_flags,flag_manual_avg,bandpass
   ```

#### Flux Density Scale Issues

**Problem**: Secondary calibrators have unreasonable flux densities.

**Solution**:

1. Check the antenna selection for fluxscale
2. Try excluding problematic antennas:

   ```bash
   # In manual_avg.flags file
   echo "mode='manual' antenna='Lo'" > manual_avg.flags
   
   # Re-run calibration
   emcp -r restore_flags,flag_manual_avg,fluxscale
   ```

### Data Issues

#### Missing Fields

**Problem**: Some fields in your dataset don't appear in the pipeline output.

**Solution**: Make sure all fields are correctly specified in the inputs.ini file.

#### Bad Antennas

**Problem**: A particular antenna shows bad data throughout the observation.

**Solution**: Flag the problematic antenna:

```bash
# In manual_avg.flags file
echo "mode='manual' antenna='Mk2'" > manual_avg.flags

# Re-run with flagging
emcp -r restore_flags,flag_manual_avg,applycal_all
```

## Improving Performance

### Parallelization

Use MPI to process data in parallel:

```bash
mpicasa -n 8 python /path/to/eMERLIN_CASA_pipeline/eMERLIN_CASA_pipeline.py
```

### Optimization Tips

1. Process smaller datasets by applying channel averaging
2. Use Docker for consistent environment setup
3. Run on systems with sufficient RAM (16GB minimum, 32GB+ recommended)

## Getting Help

If you encounter issues not covered here:

1. Check the [GitHub Issues](https://github.com/e-merlin/eMERLIN_CASA_pipeline/issues) for similar problems
2. Contact e-MERLIN staff for assistance with observatory data
3. Open a new issue with:
   - Pipeline version
   - Error messages
   - Description of the dataset
   - Steps to reproduce
