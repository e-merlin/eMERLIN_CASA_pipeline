# Customizing the Pipeline

This guide explains how to modify the eMCP code to meet your specific needs and how to run your customized version.

## Understanding the Code Structure

The pipeline code is organized as follows:

- `src/eMCP/` - Main package directory containing all pipeline modules
  - `functions/` - Core processing functions
  - `utils/` - Utility functions for file handling, etc.
  - `plots/` - Visualization and plotting functions
  - `cli.py` - Command-line interface module

## Common Customization Scenarios

### 1. Modifying Pipeline Steps

If you need to adjust how a specific pipeline step works:

1. Locate the relevant function in `src/eMCP/functions/eMCP_functions.py`
2. Make your changes to the function
3. Run the pipeline using your local modified version

**Example**: Modifying the bandpass calibration parameters

```bash
# 1. First, locate the bandpass calibration function
cd /path/to/eMERLIN_CASA_pipeline
grep -r "bandpass" src/eMCP/functions/

# 2. Edit the file with your modifications
nano src/eMCP/functions/eMCP_functions.py

# Find the initial_bp_cal function and modify parameters, for example:
# Change:
#   solint = str(solint_bp)
# To:
#   solint = 'inf'
```

### 2. Adding New Functionality

To add new capabilities to the pipeline:

1. Add your new function to the appropriate module
2. Update the CLI module if you want to make it accessible from the command line

**Example**: Adding a new flagging function

```python
# In src/eMCP/functions/eMCP_functions.py

def my_custom_flagging(eMCP):
    """
    Custom flagging function for specific issues in my dataset
    """
    msfile = eMCP['msinfo']['msfile']
    print('Applying custom flagging to {}'.format(msfile))
    
    # Your custom flagging logic here
    flagdata(vis=msfile, 
             mode='manual',
             antenna='Mk2',
             timerange='2023/01/01/12:00:00~12:30:00')
    
    return eMCP
```

### 3. Running Your Modified Version

After making changes, you need to run your local version rather than the installed version:

#### Option 1: Development Installation

Install your modified version in development mode:

```bash
cd /path/to/eMERLIN_CASA_pipeline
pip install -e .
```

This creates a link to your local code, so any changes you make will be immediately available when you run `emcp`.

#### Option 2: Direct Execution

Run the pipeline directly from your modified code:

```bash
cd /path/to/eMERLIN_CASA_pipeline
python eMERLIN_CASA_pipeline.py -r calibration
```

#### Option 3: Import in Custom Script

Create your own script that imports and uses your modified functions:

```python
#!/usr/bin/env python
# my_custom_pipeline.py

import os
import sys
from eMCP.functions import eMCP_functions as em
from eMCP.utils import eMCP_utils as emutils

# Initialize pipeline
calib_dir, info_dir = emutils.create_dir_structure()
eMCP = emutils.start_eMCP_dict(info_dir)
eMCP['inputs'] = emutils.read_inputs('./inputs.ini')

# Get MS info
msfile = eMCP['inputs']['inbase'] + '.ms'
eMCP, msinfo, msfile = em.get_msinfo(eMCP, msfile)

# Run your custom function
em.my_custom_flagging(eMCP)  # Your custom function

# Continue with standard pipeline steps
eMCP = em.initial_bp_cal(eMCP, caltables)
```

Run your custom script:

```bash
python my_custom_pipeline.py
```

## Best Practices for Code Modifications

1. **Make a backup** of the original code before making changes
2. **Document your changes** with comments explaining why modifications were made
3. **Test thoroughly** before using on important data
4. **Consider contributing** useful modifications back to the main project

If you find yourself frequently making the same modifications, consider creating a configuration option in `default_params.json` rather than modifying the code directly.
