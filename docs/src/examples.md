# Examples

This page provides practical examples for different e-MERLIN CASA Pipeline usage scenarios.

## Basic Workflow Examples

### Example 1: Processing Raw Observatory Data

This example shows how to process raw FITS-IDI data from the observatory.

1. Create your working directory:

   ```bash
   mkdir my_eMERLIN_project
   cd my_eMERLIN_project
   ```

2. Copy the inputs template:

   ```bash
   cp /path/to/eMERLIN_CASA_pipeline/inputs.ini .
   ```

3. Edit the inputs.ini file with your project details:

   ```ini
   [inputs]
   inbase = my_project
   fits_path = /path/to/fits/files/
   targets = M82
   phscals = J0955+6903
   fluxcal = 1331+305
   bpcal = 1407+284
   refant = Mk2
   ```

4. Run the full pipeline:

   ```bash
   emcp -r all
   ```

### Example 2: Calibration-Only Workflow

This example shows how to run calibration steps on already imported data.

1. Ensure you have an MS file (e.g., `my_project.ms` or `my_project_avg.ms`)
2. Run the calibration steps:

   ```bash
   emcp -r calibration
   ```

### Example 3: Custom Pipeline Flow

This example shows how to run specific steps of the pipeline.

1. Run only preprocessing steps:

   ```bash
   emcp -r pre_processing
   ```

2. Run selected calibration steps:

   ```bash
   emcp -r restore_flags,init_models,bandpass,initial_gaincal,fluxscale
   ```

3. Skip certain steps during a full run:

   ```bash
   emcp -r all -s flag_aoflagger,plot_data
   ```

## Advanced Examples

### Processing Multiple Spectral Windows

When dealing with mixed-mode observations that include both continuum and spectral line data:

1. Configure spectral window handling in `default_params.json`:

   ```json
   {
     "averaging": {
       "chanavg_continuum": 4,
       "chanavg_line": 1,
       "continuum_spws": "0,1,2,3",
       "line_spws": "4,5"
     }
   }
   ```

2. Run the pipeline:

   ```bash
   emcp -r all
   ```

### Running with Parallelization

To use multiple CPU cores with MPI:

```bash
mpicasa -n 8 python /path/to/eMERLIN_CASA_pipeline/eMERLIN_CASA_pipeline.py
```

### Using Docker

Example of using the Docker container:

```bash
# Mount your data directory and run the pipeline
docker run -it --rm -v /path/to/your/data:/data emerlin/emcp:latest emcp -r all

# Run just calibration steps
docker run -it --rm -v /path/to/your/data:/data emerlin/emcp:latest emcp -r calibration
```

## Troubleshooting Examples

### Dealing with Bad Antennas

If you identify an antenna with problems:

1. Create a `manual_avg.flags` file:

   ```bash
mode='manual' antenna='Mk2'
```

2. Run with flag restoration and manual flagging:

   ```bash
   emcp -r restore_flags,flag_manual_avg,bandpass,initial_gaincal
   ```

### Re-running Calibration with Modified Parameters

To re-run calibration with modified parameters:

1. Edit `default_params.json` to change calibration parameters
2. Restore flags and run calibration again:

   ```bash
   emcp -r restore_flags,bandpass,initial_gaincal,fluxscale,bandpass_final,gaincal_final,applycal_all
   ```
