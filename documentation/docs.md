<!---
I use grip to convert to html using: grip docs.md --export --title "e-MERLIN CASA pipeline"
It can also be converted to pdf using: http://www.markdowntopdf.com/
-->

# e-MERLIN CASA pipeline
### Documentation for v1.0.0

---
# Table of contents
- [1. Minimal execution of the pipeline](#1-minimal-execution-of-the-pipeline)
- [2. How to reduce e-MERLIN data](#2-how-to-reduce-e-merlin-data)
- [3. Inputs](#3-inputs)
- [4. Procedures](#4-procedures)
     - [4.1 Pre-processing data](#41-pre-processing)
     - [4.2 Calibration](#42-calibration)
- [5. Support variables](#5-support-functions-and-variables)

<!---
- [6. Quick summary](#6-quick-summary)
--->

---

# 1. Minimal execution of the pipeline

Go to your working directory:

`cd /here/I/work/`

Download the pipeline from github [e-MERLIN CASA Pipeline](https://github.com/e-merlin/eMERLIN_CASA_pipeline):

`git clone https://github.com/e-merlin/eMERLIN_CASA_pipeline.git`

Copy the inputs.txt file from the pipeline directory to the current location:

`cp eMERLIN_CASA_pipeline/inputs.txt .`

Edit the inputs file with the location of your FITS-IDI files and set a name for the project and the list of sources to process.

You are ready to execute the pipeline. If `casa` points to casa 5.4 version:

`casa -c eMERLIN_CASA_pipeline/eMERLIN_CASA_pipeline.py -i inputs.txt`

---

# 2. How to reduce e-MERLIN data

### Data preparation
- Download the pipeline (it can be in any location).
- Download the fits-IDI files from the observatory (they can be in any unique location).
- Create the working path where you will work, and copy the `inputs.txt` to your working path.
- Edit the inputs.txt file: fill the `fits_path` (where the fits files are located) and the `inbase` (any name you want to give to your project).
- Fill in the inputs.txt file the names of your fields: `targets`, `phscals`, `fluxcal`,`bpcal`,`ptcal`.
- Leave all other parameters as default.
- A MS `<inbase>.ms` will be produced
- Data will be converted to MS and prepared. Open in a web browser the file `./weblog/index.html`.
- The weblog and the plots will be updated every time a step is finished.
- If you have a list of manual flags to apply, write them in `inputfg.flags` to flag data before it is averaged in the next step.
- If `average = 1` a new MS will be produced with name `<inbase>_avg.ms`.

### Data calibration
- You may want to include additional manual flags in `./inputfg_avg.flags`
- Run all the calibration steps by setting them to 1. You may prefer to run each of them one by one and check the output plots.
- It is a good practice to redo the calibration once you are happy with your flags and you are sure of all the steps.
- The whole calibration process can be repeat from scratch by setting all the steps in the Calibration section to 1, including `restore_flags` which will restore the flag status when the data was averaged.


---

# 3. Inputs

In the inputs.txt file there are two types of inputs: the **user inputs** and the **process inputs**. The user inputs expect a string with information about the project, like a name, source names, or location of the fits-IDI files. The process inputs expect an integer that can be 0 to not run a step, or 1 to run it.

### User inputs

```
fits_path  [str]
```
Path to the location of the fits files. Can be an absolute or relative path. A single directory must be specified. All fits and FITS files from that directory are considered.


```
inbase     [str]
```
Project name. It will be used to give names to the MS, plots and tables.

```
targets    [str]
```
Names of sources as they appear in the MS to be used as targets. Can be a comma-separated string. If more than one target is selected, the corresponding phase calibrator has to keep the same order.


```
phscals    [str]
```
Names of sources as they appear in the MS to be used as calibrators. Can be a comma-separated string. If more than one phase calibrator is selected, the order needs to match the targets. If the same phase calibrator is used for different targets, just repeat the name as many times as needed.

```
fluxcal    [str]
```
Name of source as it appears in the MS to be used as flux calibrator. Only one source accepted. 1331+305 expected. If a different source is selected you will need to set its flux density and spectral index in the `default_params.json` file.

```
bpcal      [str]
```
Names of sources as they appear in the MS to be used as calibrators. Can be a comma-separated string.

```
ptcal      [str]
```
Point-like calibrator. Names of sources as they appear in the MS to be included in calibration steps, but not used to calibrate other sources by now. For the moment, consider it for check sources.

```
refant     [str]
```
Antenna name to be used as reference antenna. Can accept a comma-separated string with a list of antennas in order of priority. If empty, the pipeline will try to search for the most suitable antennas.

**Note on the format**: For `[str]` inputs, the use of single quotes is not required. Most inputs accept a list of values in the format of a comma-separated string. Examples: `target = '1111+2222'` or `target = 1111+2222`. For multiple inputs: `phscals = '1111+2222,3333+4444,5555+6666'` is accepted.

### Process inputs

This inputs just select which of the processing steps will be executed:

- 0 means do not execute the step

- 1 means execute the step

Additionally, steps that produce calibration tables can be set to also apply the calibration up to that point, so also modifying the corrected column. A value of 2 means run the step and apply the calibration up to that step. The standard procedure is to apply everything at the end of the calibration with the step `applycal_all = 1`, so setting 2 to any previous step is useful just to check the calibration up to each step.


---

# 4. Procedures


## 4.1 Pre-processing

### 4.1.1 run_importfitsidi
Merge fits-IDI files in `fits_path` to form an MS named `inbase.ms`

First, `inbase.ms` and any other temporary file is removed if present (`inbase` is the name of your project, which is set in inputs file inputs.txt). All files ending in .fits or .FITS in the folder `fits_path` are considered for importing and will be merged. The logger shows which files have been found. Several CASA tasks are executed: `importfitsidi` is executed with `constobsid=True, scanreindexgap_s=15.0`.

Then `mstransform` will be executed with different purposes: remove the auto-correlations, produce MMS if requested, apply any a-priori averaging to the data, and separate narrow band spectral line data if present. The pipeline checks if the observations were observed in mixed mode. If so, the spectral spws will be split so `inbase.ms` will contain only the continuum data (broadband, low spectral resolution). The high resolution spw will be splitted into `inbase_sp.ms`

Finally `fixvis` is run with option `reuse=False` to fix any possible mismatch with the UVW values in the visibilities.



---
### 4.1.2 flag_aoflagger
Uses aoflagger to autoflag data using predefined strategies


It runs [aoflagger](https://sourceforge.net/p/aoflagger/wiki/Home/) on all the fields in the MS. aoflagger v2.9+ required. The aoflagger strategies are selected by field following this criteria: (1) First, check if the user has a new strategy for this field in the local folder ./aoflagger_strategies/user/, the strategy needs to match: `<field>.rfis`. (2) If not, check if the user has produced a new strategy for this field in the pipeline folder. So will search for `pipeline_path+'aoflagger_strategies/default/<field>.rfis`, (3) If nothing is found, just use the default strategy. For example for field 1234+5678 the prioritization for searching strategies is:

 1. ./aoflagger_strategies/user/1234+5678.rfis
 2. /path/to/pipeline/aoflagger_strategies/default/1234+5678.rfis
 3. /path/to/pipeline/aoflagger_strategies/default/default_faint.rfis

By default this task will flag all the spw for each field in one go. That means that it will try to read data from all the spw for each field at the same time. If your data set is too large (compared to the available memory) it is recommended to change the default option so each spw is load one by one, with the overhead of restarting aoflagger once per spw (and per source).


---
### 4.1.3 flag_apriori
Applies a-priori standard flags.

Different flags are applied to the data:

 - Lo&Mk2 for all sources.
 - Edge channels for all sources, all spw, defined as spw='*:0~(nchan/128-1);(nchan-nchan/128)~(nchan-1)'.
 - Edge of the whole band: first 5% of the channels of the first spw, and last 5% of the channels of the last spw.
 - Quack 4 seconds on all sources in the MS.
 - It will try to find observatory flags (to be writen in `inputsfg.flags`). If not possible will produce the following ad-hoc quack commands:
 - Quack for 5 minutes for 1331+305, 1407+284, 0319+415 if they are present in the MS (no need to include in inputs file.
 - Quack n seconds for all targets and phasecals. The number of seconds depend on the target-phasecal separation in degrees. 20s for separation < 1deg, 25s for 1.0 <= separation < 2.0, 30s for 2.0 <= separation < 3.5, 35s for separation >= 3.5.
 

---
### 4.1.4 flag_manual
Applies flags from an external file with a list of flag commands the unaveraged dataset. It needs file `./inputfg.flags` to be located in the current directory.

Simply runs `flagdata(vis='inbase.ms', mode='list', inpfile="inputfg.flags")`. This is run on the original, unaveraged, dataset `inbase.ms`. To apply manual flags to an averaged dataset see `flag_2b_manual` below.

 For more information on this flagdata mode see [flagdata](https://casa.nrao.edu/casadocs/casa-5.1.1/global-task-list/task_flagdata/about). For details on the format of the file see [examples](https://casa.nrao.edu/casadocs/casa-5.1.1/global-task-list/task_flagdata/examples). Note from the [CASA documentation](https://casa.nrao.edu/casadocs/casa-5.1.1/global-task-list/task_flagdata/about) on the format of the file: There should be no whitespace between KEY=VALUE since the parser first breaks command lines on whitespace, then on "=". Use only one whitespace to separate the parameters (no commas). For example:

```
mode='manual' field='1331+305' antenna='' timerange='10:00:00~10:11:30'
mode='manual' field='' antenna='' timerange='' spw='0:0~30'
mode='manual' field='' antenna='Mk2' timerange='09:05:00~16:27:00'
mode='manual' field='1258-2219' antenna='' timerange='12:57:01~12:59:59'
mode='quack' field='1258-2219,1309-2322' quackinterval=24.
```

---
### 4.1.5 average
Split dataset and average to reduce data volume.

It will create a new MS `inbase_avg.ms` (or inbase_avg.mms if working with mms files). It will remove previous inbase_avg.ms file if it exists. Only the fields selected in the inputs file will be included in inbase_avg.ms. 

**Important** when the inbase_avg.ms file is produced, all the steps of the pipeline after this one will always work on inbase_avg.ms only, and not the original dataset. If the averaged dataset does not exist, the pipeline will always use the unaveraged dataset `inbase.ms`.

#### Phase-shift option

There is an option (to be selected from the `default_params.json` file) to shift the produce new pointings for the target fields before averaging. For that the pipeline needs a text file `shift_phasecenter.txt` containing one shift per line, with three comma-separated values each: field, new field name, and new position coordinates. For example:

```
0336+3218, 0336+3218_pos1, J2000 03h36m30.10s +32d18m30.0s
0336+3218, 0336+3218_pos2, J2000 03h36m30.10s +32d18m40.0s
0332+3205, 0336+3218_pos1, J2000 03h32m28.30s +32d05m46.0s
```

Spaces will be ignored except inside the coordinate string. That particular example file will produce two shifts of field 0336+3218 and one shift of field 0332+3205 to the corresponding positions indicated. For each shift, the field will be splitted in a temporary MS on which the shift will be computed. Then that temporary MS will be merged to the main data set. 

WARNING1: this step will add new fields to the MS (using concat) without removing any previous field: this adds phase centers, does not update the original ones. That means that the resulting MS will have two or more **new** fields with the same time stamps and same scan number but with different phase centers. This task is not intended to correct the phase center of a field, as it will keep the original one. 

WARNING2: note that because new fields are added to the MS, if you want to use them you need to specify them in the inputs file.


---
### 4.1.6 plot_data
Produce plots of amp/phase vs time/freq for each baseline	plotms

Produces plots in png format for the visibilities in `inbase_avg.ms`. It iterates through all sources specified by the user in the inputs file, and produces plots for each baseline. All plots are stored in the output directory. The plots have some averaging:

 - Amp/Phase vs Time: time averaging is 4s, all channels averaged in each spw, colorized by spw.
 - Amp/Phase vs Frequency: time average is 300s, frequency average is 4 channels. ccolorized by correlation.


### 4.1.7 save_flags
Saves the current status of the flagging to the flag table with versionname='initialize_flags'.

Makes a copy of the current flags to the flag table with versionname `initialize_flags` using flagmanager. It will overwrite previous versions of that table. This table means to contain flags produced by aoflagger, a-priori and manual_a flags, but only if the user has saved this table after producing those flags (so it is the responsability of the user to keep track of what was applied when this step is executed).


---

## 4.2. Calibration

You can select which steps to run in the inputs file by setting the corresponding number to:

 - 0: don't run the step
 - 1: run the step and produce calibration tables and plots
 - 2: run the step and produce calibration tables and plots and apply the calibration to the data. *don't use this option if you don't understand it*.

During calibration, solutions are found using only the inner ~90% channels of each spw. So `innerchan=0.1*(nchan-nchan/512.), 0.9*(nchan-nchan/512.)`.

**Important.** Below we list the parameters used to create each table in each calibration step. `previous_cal` is a list of names of previous calibration tables that will be used to fill the `gaintable`, `gainfield`, `interp`, `spwmap`. Those parameters are defined when the table is created, and they will be used only to generate other tables when the corresponding table name is listed in `previous_cal`. To know the values of those parameters you have to check them in the description of the corresponding table.

The step applycal_all will apply all calibration: it will write `CORRECTED_DATA` column in the MS using the last caltables. In particular it will use `delay.K1,bpcal_sp.B1,allcal_p.G0,allcal_p_jitter.G0,allcal_ap.G3` for the calibrators and `delay.K1,bpcal_sp.B1,phscal_p_scan.G2,allcal_ap_scan.G3` for the targets. If you want to apply the calibration up to an intermediate step (for instance you want to apply the delay and BP caltables only to check their effects before proceeding), you can set in the inputs file a value of `2`, which means apply all previous calibration up to that step. Below you can find the **applycal** description at the end of each step that shows which tables will be applied to calibrators and targets respectively if `2` is selected for that step in the inputs file.

`calsources` is a comma-separated list of all calibrator sources.


### 4.2.1 restore_flags
Restore the status of the flagging that was saved in the flag table with versionname='initialize_flags'.

Output:
```
initialize_flags        [flag table]
```

Restores the flags saved in the flag table with versionname `initialize_flags` using flagmanager. The idea is to use this saving point to restart calibration from scratch but keeping the flags produced during the pre-process stage. Any other tables associated with the dataset will not be removed, so use flagmanager and remove them manually if needed.



---
### 4.2.2 flag_manual_avg
Applies flags from an external file with a list of flag commands the unaveraged dataset. It needs file ./inputfg_avg.flags to be located in the current directory.

Simply runs `flagdata(vis='inbase_avg.ms', mode='list', inpfile="inputfg_avg.flags")`. This is run on the averaged dataset `inbase_avg.ms`. To apply manual flags to an averaged dataset see `flag_manual` step above.

For more information on this flagdata mode see [flagdata](https://casa.nrao.edu/casadocs/casa-5.1.1/global-task-list/task_flagdata/about). For details on the format of the file see [examples](https://casa.nrao.edu/casadocs/casa-5.1.1/global-task-list/task_flagdata/examples). It needs one flag command per line, with parameters separated by one space (don't use commas!). For example:

```
mode='manual' field='1331+305' antenna='' timerange='10:00:00~10:11:30'
mode='manual' field='' antenna='' timerange='' spw='0:0~30'
mode='manual' field='' antenna='Mk2' timerange='09:05:00~16:27:00'
mode='manual' field='1258-2219' antenna='' timerange='12:57:01~12:59:59'
mode='quack' field='1258-2219,1309-2322' quackinterval=24.
```
---
### 4.2.3 init_models

The pipeline will initialize the model column for all sources different from 1331+305 using CASA `delmod`, so setting all amplitudes to 1 and all phases to 0 in the model column. For 1331+3030 it will check if the data is L band or C band, and use `setjy` to introduce the correct model of 1331+3030 (3C286) into the data column. The models can be found in `pipeline_path+'calibrator_models/'`. Only C and L band models are currently available.

---
### 4.2.5 bandpass

It runs a delay, phase, and a&p calibration before finding the combined BP table for all sources listed in `bpcals`. Only the bandpass table is used in the following steps.

Tables produced:

- bpcal_d.K0 - delay solutions
- bpcal_p.G0 - phase-only solutions
- bpcal_ap.G0 - amplitude & phase solutions
- bpcal.BP0 - initial bandpass table

The pipeline will do a first pass to create the tables. Then will apply all of them to the `bpcal` sources and run an automatic flagging step (tfcrop by default). Once the data are cleaned from bad visibilities, it will do a second pass to find the four tables again.


---
### 4.2.6 initial_gaincal

It will use table bpcal.BP0 and solve for delays, phases and amplitude&phases on all calibrators.

Tables produced:

- allcal_d.K1 - delay solutions
- allcal_p.G1 - phase-only solutions
- allcal_ap.G1 - amplitude & phase solutions

The pipeline will do a first pass to create the tables. Then will apply all of them to the calibrator sources and run an automatic flagging step (tfcrop by default). Once the data are cleaned from bad visibilities, it will do a second pass to find the three tables again.

---
### 4.2.7 fluxscale

Runs CASA `fluxscale` to bootstrap the flux density scale of 1331+3030 using its model, and forward the corrections to all other sources, updating their model column. Also, a corrected `allcal_ap.G1_fluxscale` table is derived from the previous amplitude calibration.

CASA `fluxscale`

|Parameter      | Value                                     |
|-------------- | ----------------------------------------- |
| reference     | fluxcal                                   |
| transfer      | calibrators (except fluxcal)              |
| caltable      | allcal_ap.G1                              |
| fluxtable     | allcal_ap.G1_fluxscaled                   |
| listfile      | ./calib/allcal_ap.G1_fluxscaled_fluxes.txt|
|               |                                           |


Then, the pipeline will run the script `dfluxpy` to find the correction factor eMfactor. This is a correction factor for the flux density of 1331+3030 (3C286) needed because the source is slightly resolved by the shortest baseline of e-MERLIN. The task will check the shortest baseline and the observation frequency, and scale the results accordingly. The values reported by the logger in eMCP.log are already corrected by this factor. The file allcal_ap.G1_fluxscaled_fluxes.txt is not corrected by this factor, but a warning note is included in the file.

The step can accept an external model image for any source (typically the phase calibrator). For each source `fieldname` it will search for `./source_models/<fieldname>.model.tt0` and `./source_models/<fieldname>.model.tt1`. If both are present, the models will be scaled by eMcalflux/flux_in_model. eMcalflux is the flux derived by task `fluxscale` (see above), and flux_in_model is the sum of pixels in the model.tt0 image. Both the tt0 and tt1 model images will be scaled by that factor and the new models specific for the current observation will be placed in ./source_models. If no models are found for a source, a point-like model with the flux density and spectral index derived by `fluxscale` will be used. Finally, the model column in the MS will be updated using `ft` if there model images (tt0 and tt1) are available in `./source/models/`, of `setjy` if there are no models and a point-like model is assumed. In both cases the model column should contain spectral index information.


---
### 4.2.8 bandpass_final

Tables produced:

- bpcal.BP2 - final bandpass correction (includes spectral information)

Recalculate the BP table. Now the model columns contain the actual flux density of the calibrators, so the BP table will include spectral index information.


---
### 4.2.9 gaincal_final

With updated models on all calibrators and the final bandpass (BP2) we can recompute phase and amplitude solutions. Note delays from `initial_gaincal` are used.

Tables produced:

- allcal_p.G3
- allcal_ap.G3
- phscal_p_scan.G3
- phscal_ap_scan.G3

The first two use a short solution interval (by default `int` for phases and `32s` for a&p) and the later two produce scan-averaged solutions to correct the targets. The scan-averaged solutions are only found for phase-reference sources, but not for the bandpass or flux calibrator.

---
### 4.2.10 applycal_all
Two runs are executed. One to correct calibrators, using their own solutions when relevant. A second correction is executed for each target, using the solutions from the corresponding phase reference calibrator.


---
### 4.2.11 flag_target

After calibration, we can run `flagdata` in tfcrop mode to remove RFI from the target fields.

---
### 4.2.13 plot_corrected

Produce plots of amp/phase vs time/freq for each baseline with plotms using corrected data column. Also produce Amp/phase vs UVwave plots of corrected visibilities and model.


Produces plots in png format for the visibilities in `inbase_avg.ms`. It iterates through all sources specified by the user in the inputs file, and produces plots for each baseline. All plots are stored in the output directory. The plots have some averaging:
 
 - Amp/Phase vs Time: time averaging is 4s, all channels averaged in each spw, colorized by spw.
 - Amp/Phase vs Frequency: time average is 300s, frequency average is 4 channels. ccolorized by correlation.
 - Amp/Phase vs UVwave: time average is 600s, channel average is 1/16 of total channels


---

# 5. Support functions and variables

### msinfo [dict]
The pipeline stores all available information on the dataset and the steps executed and all the parameters of each step in a dictionary. The dictionary is saved in pickle in file in `./weblog/info/eMCP_info.pkl`. It can be read with picke using:

> `msinfo = pickle.load(open('weblog/info/eMCP_info.pkl', 'rb'))`

The contents are also in human-readable txt file: `weblog/info/eMCP_info.txt`, which can also be accessed from the weblog, tab `Pipeline info`.

Example:

```
  pipeline_path     : /mirror2/scratch/jmoldon/test_pipeline/test_complete5b/eMERLIN_CASA_pipeline/
  casa_version      : 5.4.0
  pipeline_version  : v0.10.21
  msfile            : DD6001_C_001_20171220.ms
  is_mixed_mode     : False
  img_stats
    1109-1235       : [0.0066674887202680111, 0.00013494912213562872, 0.0]
    1118-1232       : [0.710640549659729, 0.0050747106203501112, -0.94189724325111968]
    1107-1226       : [0.056878048926591873, 0.0011088359273201554, 0.0]
  inputs
    fits_path       : ../data/
    inbase          : DD6001_C_001_20171220
    targets         : 1109-1235,1107-1226
    phscals         : 1118-1232,1118-1232
    fluxcal         : 1331+305
    bpcal           : 1407+284
    ptcal           : 1407+284
    refant          : 
    run_importfits  : 0
    flag_aoflagger  : 0
    flag_apriori    : 0
    flag_manual     : 0
    average         : 0
    plot_data       : 0
    save_flags      : 0
    restore_flags   : 0
    flag_manual_avg : 0
    init_models     : 0
    bandpass        : 0
    initial_gaincal : 0
    fluxscale       : 0
    bandpass_final  : 0
    gaincal_final   : 0
    applycal_all    : 0
    flag_target     : 0
    plot_corrected  : 0
    first_images    : 0
  defaults
    aoflagger
        separate_bands: False
        fields      : all
        run         : auto
    bp_apply_mid
        apply_targets: []
        apply_calibrators: ['bpcal_d.K0', 'bpcal_p.G0', 'bpcal_ap.G0', 'bpcal.BP0']
    flag_manual_avg
        Lo_datacolumn: data
        Lo_threshold: 0.5
        Lo_min_scans: 
        Lo_dropout  : 
        Lo_spws     : ['3']
        Lo_useflags : True
    plot_data
        num_proc    : 1
    applycal_all
        apply_narrow_targets: ['allcal_d.K1', 'narrow_bpcal.BP2', 'phscal_p_scan.G3', 'phscal_ap_scan.G3', 'narrow_p_offset.G3']
        apply_targets: ['allcal_d.K1', 'bpcal.BP2', 'phscal_p_scan.G3', 'phscal_ap_scan.G3']
        apply_narrow_calibrators: ['allcal_d.K1', 'narrow_bpcal.BP2', 'allcal_p.G3', 'allcal_ap.G3', 'narrow_p_offset.G3']
        apply_calibrators: ['allcal_d.K1', 'bpcal.BP2', 'allcal_p.G3', 'allcal_ap.G3']
        statwt_timebin: 0.001s
        run_statwt  : True
    average
        timebin     : 4s
        antenna     : 
        scan        : 
        timerange   : 
        datacolumn  : data
        field       : 
        chanbin     : 4
        shift_phasecenter: False
...
...

```
### caltables [dict]

This is a dictionary that is stored in the file `./calib/caltables.pkl`. It contains information on all the calibration tables produced by the pipeline. The pipeline saves cumulative pickle files after each step as ./calib/caltables_<step>.pkl. This is just for security and it is not used by the pipeline. The caltables variable can be load with python using:

> `caltables = pickle.load(open('./calib/caltables.pkl', 'rb'))`


Example:

```
'bpcal_d.K0': {'calmode': 'p',
                'combine': 'spw',
                'field': '1407+284,1258-2219',
                'gainfield': '1407+284,1258-2219',
                'gaintype': 'K',
                'interp': 'linear',
                'name': 'bpcal_d.K0',
                'previous_cal': [],
                'solint': '180s',
                'spw': '*:13~115',
                'spwmap': [0, 0, 0, 0],
                'table': './calib/CY5003_24_C_20170922_bpcal_d.K0'},
 'bpcal_p.G0': {'calmode': 'p',
                'combine': '',
                'field': '1407+284,1258-2219',
                'gainfield': '1407+284,1258-2219',
                'gaintype': 'G',
                'interp': 'linear',
                'name': 'bpcal_p.G0',
                'previous_cal': ['bpcal_d.K0'],
                'solint': 'int',
                'spw': '*:13~115',
                'spwmap': [],
                'table': './calib/CY5003_24_C_20170922_bpcal_p.G0'},

... 
...
```


<!---

---

# 6. Quick summary

### Pre-processing summary
|Procedure        |Summary                                                     | CASA/external tasks                       |Inputs                   |Outputs        |Notes                                                                                          |
|-----------------|------------------------------------------------------------|-------------------------------------------|-------------------------|---------------|-----------------------------------------------------------------------------------------------|
|run_importfitsidi| Concatenate all fits-IDI files in a folder to a MS         | importfitsidi, fixvis, flagdata (autocorr)|fits_path                |inbase.ms      | Also produces inbase.ms.listobs.txt                                                               |
|summary_weblog   | Produces weblog with basic project information and plots   | get_data, plotms                          |                         |               | Info (sources, spw, channels, time, etc), elevation vs time plot and uvcov plots.             |
| hanning         | Run hanning smoothing                                      | mstransform                               |                         |               | Hanning smoothing helps with RFI removal. Recommended for L band data.                        |                                                                                                                                                                              | ms2mms          | Convert to mms format to run tasks in parallel             | mstransform                               |                         |inbase.mms     | If this step is executed, inbase.ms will be removed and all other steps will run on inbase.mms|
| flag0_aoflagger | Autoflag data using predefined strategies                  | aoflagger                                 |                         |               | Accepts user strategies per field. Can take very long if low memory available                 |
| flag1_apriori   | Standard flags: Lo-Mk2, edge channels, slewing, Lo dropouts| flagdata                                  |sources, Lo_dropout_scans|               |                                                                                               |
| flag_2a_manual  | Applies flags from external file with a list of flag commands | flagdata                               |manual_flags_a           |               | Applies flag commands prepared by the user in external file to **unaveraged** data            |
| average_1       | Split dataset and average to reduce data volume            | split                                     |sources                  |inbase_avg.ms  |                                                                                               |
| plot_data       | Produce plots of amp/phase vs time/freq for each baseline  | plotms                                    |                         |plots (several)| Plots DATA column only.                                                                       |



### Calibration summary

The calibration steps require `targets`, `phscals`, `fluxcal`, `bpcal`. All calibration steps that produce calibration tables also produce plots of those tables that are saved to `./plots/caltables/`.

|Procedure      |Summary                                                                | CASA/external tasks |Output table(s)                                                | Notes                                                                                                                        |
|---------------|---------------------------------------------------------------------- |---------------------|---------------------------------------------------------------|------------------------------------------------------------------------------- |
|flag_2b_manual | Applies flags from external file with a list of flag commands         |flagdata             |manual_flags_b                                                 |Applies flag commands prepared by the user in external file to **averaged** data|
|init_models    |Initializes MODEL column                                               |setjy                |                                                               |Initialize all other sources to amp=1, phase=0. 1331+305 uses 3C286 models at L and C bands                                   |
|bandpass_0     |Initial calibration of bandpass calibrator and creates initial BP table|gaincal, bandpass    |bpcal_d.K0, bpcal_p.G0, bpcal_ap.G1, bpcal.B0                  |Performs delay, phase and a&p calibration on bpcal sources. Runs `bandpass` with solnorm=True                                 |
|flag_3_tfcropBP|Autoflag BP-corrected data                                             |flagdata             |                                                               |Will apply BP cal if not done before in the current pipeline run. Runs fagdata with mode='TFCROP'                             |
|delay          |Baseline-based delay calibration of all calibrators                    |gaincal              |delay.K1                                                       |gaintype='k' with 300s solution interval. Can use run_fringefit when fringefit is available.                                  |
|gain_0_p_ap    |Phase and amplitude calibration against time for all calibrators       |gaincal              |allcal_p.G0, allcal_p_jitter.G0, allcal_ap.G1, phscal_p_scan.G2|Two phase calibrations, a&p. Also phase calibration solint='inf' on phscals for phase referencing                             |
|fluxscale      |Find abolute flux density scale by bootstraping to 1331+305 (3C286)    |fluxscale            |allcal_ap.G1_fluxscaled, fluxes.txt                            |Will try to not include Lo, and then De, while keeping at least 4 station
|bandpass_1_sp  |Re-calculate BP table including real flux densities and sp index       |bandpass             |bpcal_sp.B1                                                    |Same phase calibration tables and allcal_ap.G1_fluxscaled with solnorm=False                                                  |
|gain_1_amp_sp  |Re-calculate amplitude calib including spectral index                  |gaincal              |allcal_ap.G3, allcal_ap_scan.G3                                |Uses original phase calibration and new bandpass to find a&p solutions. Also a&p on phscals solint='inf' for phase referencing|
|applycal_all   |Apply calibration and fill corrected_data column                       |applycal             |                                                               |Apply calibration to all sources. Targets get scan-averaged solutions for p and a&p.                                          |
|flag_4_rflag   |Autoflag data based on corrected_data column                           |flagdata             |                                                               |Runs fagdata with mode='RFLAG'                                                                                                |
|plot_corrected |Produce plots of amp/phase vs time/freq for each baseline              |               plotms|                                                               |Plots CORRECTED data column                                                                                                   |
|weblog         |Produces weblog with project information and plots                     |                     |                                                               |Now includes uncalibrated/calibrated visibilities, tables, and uvplots                                                        |

-->
