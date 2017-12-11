<!---
I use grip to convert to html using: grip docs.md --export --title "e-MERLIN CASA pipeline"
It can also be converted to pdf using: http://www.markdowntopdf.com/
-->

# e-MERLIN CASA pipeline
### Documentation for v0.6

---
# Table of contents
- [1. How to run the pipeline](#1-how-to-run-the-pipeline)
- [2. How to reduce e-MERLIN data](#2-how-to-reduce-e-merlin-data)
- [3. Inputs](#3-inputs)
- [4. Procedures](#4-procedures)
     - [4.1 Pre-processing data](#41-pre-processing)
     - [4.2 Calibration](#42-calibration)
- [5. Support functions](#5-support-functions-and-variables)
- [6. Quick summary](#6-quick-summary)


---

# 1. How to run the pipeline


Download the pipeline from github [e-MERLIN CASA Pipeline](https://github.com/e-merlin/CASA_eMERLIN_pipeline).

To run the pipeline simply do:
```
casa -c /path/to/pipeline/eMERLIN_CASA_pipeline.py -i <input file>
```

To run the parallelized version using MPI in CASA you can use:
```
mpicasa -n <num_cores> casa -c /path/to/pipeline/eMERLIN_CASA_pipeline.py -i <input file>
```

To execute the pipeline from within CASA:
~~~~
run_in_casa = True
pipeline_path = '/path/to/pipeline_path/'   # You need to define this variable explicitly
execfile(pipeline_path + 'eMERLIN_CASA_pipeline.py')
inputs, msinfo = run_pipeline(inputs_path=<input file>)
~~~~

---

# 2. How to reduce e-MERLIN data

### Data preparation
- Download the pipeline (it can be in any location).
- Download the fits-IDI files from the observatory (they can be in any unique location).
- Create the working path where you will work, and copy the `inputs.txt` to your working path.
- Edit the inputs.txt file: fill the `fits_path` (where the fits files are located) and the `inbase` (any name you want to give to your project).
- Leave all other parameters as default, and only use the steps `run_importfits=1` and `summary_weblog=1`.
- Data will be converted to MS and prepared. Open in a web browser the file `./weblog/index.html`.
- Check the listobs file in the 'Observation summary' tab and fill in the inputs.txt file the fields `targets`, `phscals`, `fluxcal`,`bpcal`,`ptcal`.
- Select a reference antenna `refant`. If you are not sure set `refant=''` and rerun the `summary_weblog` step alone. The pipeline will try to suggest a list of the best reference antennas to use.
- Run the rest of the pre-processing steps depending on your needs. If you have a list of manual flags to apply, set `flag_2a_manual=1` and remember to specify where that file is by setting `manual_flags_a` in the user inputs.
- In most cases it is recommended to run `average_1=1` to split the date to a new averaged dataset `inbase_avg.ms`.

### Data calibration
- Prepare the file `manual_flags_b` with flags commands based on the plots produced in the previous section (or your own data exploration). The external file should be set in `manual_flags_b` in the user inputs.
- Run all the calibration steps by setting them to 1. You may prefer to run each of them one by one and check the output plots
- You can always select the step `weblog=1`. It will just update the weblog with any new plots available.
- Improve the `manual_flags_a` file to have more detailed flagging as you proceed with the calibration.
- It is a good practice to redo the calibration once you are happy with your flags and you are sure of all the steps. For that, you can rerun `average_1` to produce from scratch the `inbase_avg.ms` dataset and repeat all the calibration steps.



---

# 3. Inputs

There are two types of inputs: the **user inputs** and the **process inputs**. The user inputs expect a string that will depend on the project, sources, antennas and external files. The process inputs expect an integer that can be 0 to not run a step, or 1 to run it. Additionally, steps that produce calibration tables can be set to also apply the calibration to the data modifying the corrected column. So for `bandpass_0`, `delay`, `gain_0_p_ap`, `fluxscale`, `bandpass_1_sp` and `gain_1_amp_sp`, a value of 2 means run the step and apply the calibration up to that step. This is useful to check the calibration up to each step, but you can also use `applycal_all=1` to apply everything when all tables are produced.


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
Name of source as it appears in the MS to be used as flux calibrator. Only one source accepted. 1331+305 expected. If a different source is selected the `fluxscale` step will not work properly.


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
Antenna name to be used as reference antenna. Can accept a comma-separated string, but CASA can only manage a single value for now (there are plans to allow `gaincal` to use a prioritized list. If empty, the pipeline will try to search for the most suitable antennas.


```
Lo_dropout_scans    [str]
```
Comma-separated list of scans in which the Lo telescope was not observing the phase calibrator(s). This has two effects: those scans will be flagged in the `flag_1_apriori` step, and also each calibration table created will be edited to remove solutions for antenna Lo for the phase calibrator scans listed in `Lo_dropout_scans`. Example: `Lo_dropout_scans = '4,8,12,16,20,24,28'`

```
manual_flags_a      [str]
```
Path to an external file that contains a list of flag commands, one per line. The flags will be applied to the unaveraged dataset `inbase.ms` when step `flag_2a_manual` is enabled. The file will be read by CASA task flagdata using `mode='list'` and `inpfile` will be set to `manual_flags_a`. Note from the [CASA documentation](https://casa.nrao.edu/casadocs/casa-5.1.1/global-task-list/task_flagdata/about) on the format of the file: There should be no whitespace between KEY=VALUE since the parser first breaks command lines on whitespace, then on "=". Use only one whitespace to separate the parameters (no commas).

```
manual_flags_b      [str]
```
Same as `manual_flags_a` but will be applied to the averaged data set `inbase_avg.ms` when step `flag_2b_manual` is enabled.


**Note on the format**: For `[str]` inputs, the use of single quotes is not required. Most inputs accept a list of values in the format of a comma-separated string. Examples: `target = '1111+2222'` or `target = 1111+2222`. For multiple inputs: `phscals = '1111+2222,3333+4444,5555+6666'` is accepted.

---

# 4. Procedures


## 4.1 Pre-processing

### 4.1.1 run_importfitsidi
Merge fits-IDI files in `fits_path` to form an MS named `inbase.ms`

Inputs parameters needed:
```
fits_path               [str]
inbase                  [str]
```

Output:
```
inbase.ms               [MS]
inbase.ms.listobs.txt   [txt]
inbase.ms.sp0           [MS] (only for mixed mode observations)
```

First, `inbase.ms` is removed if present (`inbase` is the parameter in inputs file). All files ending in .fits or .FITS in the folder `fits_path` are considered for importing and will be merged. The logger shows which files have been found. Several CASA tasks are executed: `importfitsidi` is executed with `constobsid=True, scanreindexgap_s=15.0`. The dataset can be averaged in time at this point by setting the `run_importfits = n` with n > 1. The dataset will be averaged with a timebin of `n`. 

Then `uvfix` with option `reuse=False` is run to correct UVW values and the output MS replaces `inbase.ms`. Then `flagdata` with `mode='manual',autocorr=True` to remove autocorrelations. Finally, `listobs` is run and a new listobs file is created in `inbase.ms.listobs.txt`.

At this point, the pipeline checks if the observations were observed in mixed mode. If so, the spectral spws will be split and `inbase.ms` will contain only the continuum data (broadband, low spectral resolution). The high resolution spw will be splitted in new MS ending with sp0, sp1, etc.


---
### 4.1.2 summary_weblog
Retrieve basic information from the original MS, produce some plots and compile everything in the weblog.

Inputs parameters needed:
```
inbase                  [str]
```

Output:
```
inbase.ms.msinfo.pkl    [pickle dictionary]
```
Runs `get_msinfo` (see section 3, below) to get information from the inputs file and directly from the MS. It prepares plots for elevation vs time and uvcov for each individual sources. It produces a weblog with the available information. If no `refant` is selected, or the antenna selected is not in the MS, a procedure will try to guess a good refant to use. Please update the inputs file with a suitable refant.

---
### 4.1.3 hanning
Just runs CASA hanning smoothing.

Inputs parameters needed:
```
inbase                  [str]
```

Output:
```
inbase.ms               [MS]
```

Runs mstransform with mode='hanning' on DATA column. It produces a new MS, which substitutes the `inbase.ms`. It works with MS or MMS. This step is optional but recommended for L band observations. It is probably not needed for C band datasets.



---
### 4.1.4 ms2mms
Optional step to transform original MS to an MMS.

Inputs parameters needed:
```
inbase                  [str]
```

Output:
```
inbase.mms               [MMS]
```

This allows CASA to run tasks on the MultiMeasurementSet in parallel transparently. Takes time to convert but later makes calibration and imaging faster. To be tested in detail. Uses CASA task `partition` with default parameters. **Important**: if ms2mms is executed, it will delete `inbase.ms` and all the following steps will be run on the new `inbase.mms`.


---
### 4.1.5 flag_0_aoflagger
Uses aoflagger to autoflag data using predefined strategies

Inputs parameters needed:
```
inbase                  [str]
```

It runs [aoflagger](https://sourceforge.net/p/aoflagger/wiki/Home/) on all the fields in the MS. aoflagger v2.9+ required. The pipeline might work with previous versions in some situations, but it is not recommended and probably will fail The strategies are selected by field following this criteria: (1) First, check if the user has a new strategy for this field in the local folder ./aoflagger_strategies/user/, the strategy needs to match: `<field>.rfis`. (2) If not, check if the user has produced a new strategy for this field in the pipeline folder. So will search for `pipeline_path+'aoflagger_strategies/default/<field>.rfis`, (3) If nothing is found, just use the default strategy. For example for field 1234+5678 the prioritization for searching strategies is:

 1. ./aoflagger_strategies/user/1234+5678.rfis
 2. /path/to/pipeline/aoflagger_strategies/default/1234+5678.rfis
 3. /path/to/pipeline/aoflagger_strategies/default/default_faint.rfis



---
### 4.1.6 flag_1_apriori
Applies a-priori standard flags.

Inputs parameters needed:
```
inbase                  [str]
targets                 [str]
phscals                 [str]
fluxcal                 [str]
bpcal                   [str]
ptcal                   [str]
Lo_dropout_scans        [str]
```
Different flags are applied to the data:

 - Lo&Mk2 for all sources.
 - Edge channels for all sources, all spw, defined as spw='*:0~(nchan/128-1);(nchan-nchan/128)~(nchan-1)'.
 - Quack for 5 minutes for 1331+305, 1407+284, 0319+415 if they are present.
 - Quack n seconds for all targets and phasecals. The number of seconds depend on the target-phasecal separation in degrees. 20s for separation < 1deg, 25s for 1.0 <= separation < 2.0, 30s for 2.0 <= separation < 3.5, 35s for separation >= 3.5.
 - Scans in `Lo_dropout_scans` for antenna = 'Lo' and fields in `phscals`.


---
### 4.1.7 flag_2a_manual
Applies flags from an external file with a list of flag commands the unaveraged dataset.

Inputs parameters needed:

```
inbase                  [str]
manual_flags_a          [str, path to file]
```

Simply runs `flagdata(vis='inbase.ms', mode='list', inpfile=manual_flags_a)`. This is run on the original, unaveraged, dataset `inbase.ms`. To apply manual flags to an averaged dataset see `flag_2b_manual` below.

 For more information on this flagdata mode see [flagdata](https://casa.nrao.edu/casadocs/casa-5.1.1/global-task-list/task_flagdata/about). For details on the format of the file see [examples](https://casa.nrao.edu/casadocs/casa-5.1.1/global-task-list/task_flagdata/examples). It needs one flag command per line, with parameters separated by one space (don't use commas!). For example:

```
mode='manual' field='1331+305' antenna='' timerange='10:00:00~10:11:30'
mode='manual' field='' antenna='' timerange='' spw='0:0~30'
mode='manual' field='' antenna='Mk2' timerange='09:05:00~16:27:00'
mode='manual' field='1258-2219' antenna='' timerange='12:57:01~12:59:59'
mode='quack' field='1258-2219,1309-2322' quackinterval=24.
```
---
### 4.1.8 average_1
Split dataset and average to reduce data volume.

Inputs parameters needed:
```
inbase                  [str]
targets                 [str]
phscals                 [str]
fluxcal                 [str]
bpcal                   [str]
ptcal                   [str]
```

Output:
```
inbase_avg.ms               [MS]
inbase_avg.ms.listobs.txt   [txt]
```

It will create a new MS `inbase_avg.ms` (or inbase_avg.mms if working with mms files). It will remove previous inbase_avg.ms file if it exists. Only the fields selected in the inputs file will be included in inbase_avg.ms. Data is averaged using `width=4` and `timebin='1s'` (This will be `2s` when CASA gaincal fixes its bug related to VisibilityIterator2). Optionally you can use a value different than `1` in the inputs file, and that number will be the timebin value in seconds. `keepflags=False` in this version.

---
**Important** when the inbase_avg.ms file is produced, all the steps of the pipeline after this one will always work on inbase_avg.ms only, and not the original dataset. At this point, `get_msinfo` is executed again, and the `msinfo` dictionary is created and saved in `inbase_avg.ms.msinfo.pkl`. If the averaged dataset does not exist, the pipeline will always use the unaveraged dataset `inbase.ms`.


---
### 4.1.9 plot_data
Produce plots of amp/phase vs time/freq for each baseline	plotms

Inputs parameters needed:
```
inbase                  [str]
targets                 [str]
phscals                 [str]
fluxcal                 [str]
bpcal                   [str]
ptcal                   [str]
```

Output:
```
./plots/plots_data      [directory with multiple plots in png format]
```

Produces plots in png format for the visibilities in `inbase_avg.ms`. It iterates through all sources specified by the user in the inputs file, and produces plots for each baseline. All plots are stored in the output directory. The plots have some averaging:

 - Amp/Phase vs Time: all channels averaged in each spw, colorized by spw.
 - Amp/Phase vs Frequency: averaged every 5 min, colorized by correlation.


### 4.1.10 save_flags
Saves the current status of the flagging to the flag table with versionname='initialize_flags'.

Output:
```
initialize_flags        [flag table]
```

Makes a copy of the current flags to the flag table with versionname `initialize_flags` using flagmanager. It will overwrite previous versions of that table. This table means to contain flags produced by aoflagger, a-priori and manual_a flags, but only if the user has saved this table after producing those flags (so it is the responsability of the user to keep track of what was applied when this step is executed).


---

## 4.2. Calibration

All calibration steps require the following input parameters:
```
inbase                  [str]
targets                 [str]
phscals                 [str]
fluxcal                 [str]
bpcal                   [str]
ptcal                   [str]
refant                  [str]
Lo_dropout_scans        [str]
```

You can select which steps to run in the inputs file by setting the corresponding number to:

 - 0: don't run the step
 - 1: run the step and produce calibration tables and plots
 - 2: run the step and produce calibration tables and plots and apply the calibration to the data.

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
### 4.2.2 flag_2b_manual
Applies flags from an external file with a list of flag commands to the averaged dataset.

Inputs parameters needed:

```
manual_flags_b          [str, path to file]
```

Simply runs `flagdata(vis='inbase_avg.ms', mode='list', inpfile=manual_flags_b)`. This is run on the averaged dataset `inbase_avg.ms`. To apply manual flags to an averaged dataset see `flag_2b_manual` below.

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

The pipeline will initialize the model column for all sources different from 1331+305 using CASA `delmod`, so setting all amplitudes to 1 and all phases to 0 in the model column. For 1331+305 it will check if the data is L band or C band, and use `setjy` to introduce the correct model of 1331+305 (3C286) into the data column. The models can be found in `pipeline_path+'calibrator_models/'`. Only C and L band models available.

---
### 4.2.5 bandpass_0

It runs a delay, phase, and a&p calibration before finding the combined BP table for all sources listed in `bpcals`.

**bpcal_d.K0** - Delay calibration of bandpass calibrator(s)

|Table parameter| Value                    |
|---------------| -------------------------|
| name          | bpcal_d.K0               |
| field         | bpcal                    |
| previous_cal  | []                       |
| solint        | 180s                     |
| gaintype      | K                        |
| calmode       | p                        |
| spw           | *:innerchan              |
| combine       | spw                      |
| table         | ./calib/inbase_bpcal_d.K0|
| gainfield     | bpcal                    |
| interp        | linear                   |
| spwmap        | [0]*num_spw              |



**bpcal_p.G0** - Phase calibration of bandpass calibrator(s)

|Table parameter| Value                    |
|---------------| -------------------------|
| name          | bpcal_p.G0               |
| field         | bpcal                    |
| previous_cal  | ['bpcal_d.K0']           |
| solint        | int                      |
| gaintype      | G                        |
| calmode       | p                        |
| spw           | *:innerchan              |
| combine       |                          |
| table         | ./calib/inbase_bpcal_p.G0|
| gainfield     | bpcal                    |
| interp        | linear                   |
| spwmap        | []                       |

**bpcal_ap.G1** - Amplitude and phase calibration of bandpass calibrator(s)

|Table parameter| Value                       |
|---------------| --------------------------- |
| name          | bpcal_p.G0                  |
| field         | bpcal                       |
| previous_cal  | ['bpcal_d.K0', 'bpcal_p.G0']|
| solint        | 32s                         |
| gaintype      | G                           |
| calmode       | ap                          |
| spw           | *:innerchan                 |
| combine       |                             |
| table         | ./calib/inbase_bpcal_ap.G1  |
| gainfield     | bpcal                       |
| interp        | linear                      |
| spwmap        | []                          |

**bpcal.B1** - Initial bandpass table

|Table parameter| Value                                      |
|-------------- | ------------------------------------------ |
| name          | bpcal.B0                                   |
| field         | bpcal                                      |
| previous_cal  | ['bpcal_d.K0', 'bpcal_p.G0', 'bpcal_ap.G1']|
| solint        | inf                                        |
| spw           |                                            |
| uvrange       |                                            |
| combine       | field,scan                                 |
| table         | ./calib/inbase_bpcal.B0                    |
| gainfield     | bpcal                                      |
| interp        | nearest,linear                             |
| spwmap        | []                                         |

**Applycal:**

 - On calibrators: ['bpcal.B0']
 - On targets: ['bpcal.B0']


---
### 4.2.5 flag_3_tfcropBP

This task needs the data to be bandpass corrected. So first of all, it will apply the table ['bpcal.B0'] to all sources.

CASA `flagdata`

|Parameter      | Value  |
|-------------- | ------ |
| mode          | tfcrop |
| correlation   | ABS_ALL|
| ntime         | 90min  |
| combinescans  | True   |
| datacolumn    | DATA   |
| winsize       | 3      |
| timecutoff    | 3.6    |
| freqcutoff    | 3.6    |
| maxnpieces    | 2      |
| usewindowstats| sum    |
| halfwin       | 3      |
| extendflags   | True   |
| action        | apply  |

---
### 4.2.6 delay

**bpcal_d.K0** - Delay calibration of all calibrators.

|Table parameter| Value                  |
|---------------| -----------------------|
| name          | delay.K1               |
| field         | calsources             |
| previous_cal  | ['bpcal.B0']           |
| solint        | 300s                   |
| gaintype      | K                      |
| calmode       | p                      |
| spw           | *:innerchan            |
| combine       | spw                    |
| table         | ./calib/inbase_delay.K1|
| gainfield     | calsources             |
| interp        | linear                 |
| spwmap        | [0]*num_spw            |

**Applycal:**

 - On calibrators: ['bpcal.B0','delay.K1']
 - On targets: ['bpcal.B0','delay.K1']


---
### 4.2.7 gain_0_p_ap

**allcal_p.G0** - Phase calibration of all calibrators.

|Table parameter| Value                     |
|---------------| --------------------------|
| name          | allcal_p.G0               |
| field         | calsources                |
| previous_cal  | ['delay.K1', 'bpcal.B0']  |
| solint        | 16s                       |
| gaintype      | G                         |
| calmode       | p                         |
| spw           | *:innerchan               |
| combine       |                           |
| table         | ./calib/inbase_allcal_p.G0|
| gainfield     | calsources                |
| interp        | linear                    |
| spwmap        | []                        |

**allcal_p_jitter.G0** - Short-interval phase calibration of all calibrators.

|Table parameter| Value                            |
|---------------| ---------------------------------|
| name          | allcal_p_jitter.G0               |
| field         | calsources                       |
| previous_cal  | ['delay.K1', 'bpcal.B0']         |
| solint        | 2s                               |
| gaintype      | G                                |
| calmode       | p                                |
| spw           | *:innerchan                      |
| combine       | spw                              |
| table         | ./calib/inbase_allcal_p_jitter.G0|
| gainfield     | calsources                       |
| interp        | linear                           |
| spwmap        | [0]*nchan                        |


**allcal_ap.G1** - Amplitude and phase calibration of all calibrators.

|Table parameter| Value                                                        |
|---------------| ------------------------------------------------------------ |
| name          | allcal_ap.G1                                                 |
| field         | calsources                                                   |
| previous_cal  | ['delay.K1', 'bpcal.B0', 'allcal_p.G0', 'allcal_p_jitter.G0']|
| solint        | 32s                                                          |
| gaintype      | G                                                            |
| calmode       | ap                                                           |
| spw           | *:innerchan                                                  |
| combine       |                                                              |
| table         | ./calib/inbase_allcal_ap.G1                                  |
| gainfield     | calsources                                                   |
| interp        | linear                                                       |
| spwmap        | []                                                           |


**phscal_p_scan.G2** - Scan-averaged phase calibration of all phase reference calibrators.

|Table parameter| Value                          |
|---------------| -------------------------------|
| name          | phscal_p_scan.G2               |
| field         | calsources                     |
| previous_cal  | ['delay.K1', 'bpcal.B0']       |
| solint        | inf                            |
| gaintype      | G                              |
| calmode       | p                              |
| spw           | *:innerchan                    |
| combine       |                                |
| table         | ./calib/inbase_phscal_p_scan.G2|
| gainfield     | calsources                     |
| interp        | linear                         |
| spwmap        | []                             |



**Applycal:**

 - On calibrators: ['delay.K1','allcal_p.G0', 'allcal_p_jitter.G0', 'allcal_ap.G1','bpcal.B0']
 - On targets: ['delay.K1','phscal_p_scan.G2','allcal_ap.G1','bpcal.B0']



---
### 4.2.8 fluxscale

Runs CASA `fluxscale` to bootstrap the flux density scale of 1331+305 using its model, and forward the corrections to all other sources, updating their model column. Also, a corrected _fluxscale table is derived from the previous amplitude calibration.

CASA `fluxscale`

|Parameter      | Value                                     |
|-------------- | ----------------------------------------- |
| reference     | fluxcal                                   |
| transfer      | calibrators (except fluxcal)              |
| antenna       | anten_for_flux*                           |
| caltable      | allcal_ap.G1                              |
| fluxtable     | allcal_ap.G1_fluxscaled                   |
| listfile      | ./calib/allcal_ap.G1_fluxscaled_fluxes.txt|
|               |                                           |

*anten_for_flux: selection of antennas to use for amplitude scale. Will try to remove Lo and De from the fluxscale determination, but only if there are enough antennas to have at least 4 antennas.


Then, the pipeline will run the script `dfluxpy` to find the correction factor eMfactor. This is a correction factor for the flux density of 1331+305 (3C286) needed because the source is slightly resolved by the shortest baseline of e-MERLIN. The task will check the shortest baseline and the observation frequency, and scale the results accordingly. The values reported by the logger in eMCP.log are already corrected by this factor. The file allcal_ap.G1_fluxscaled_fluxes.txt is not corrected by this factor, but a warning note is included in the file.

Finally, CASA `setjy` is run for each calibrator source (except the fluxcal) to include the new flux density into the model column for each spw.


**Applycal:**

 - On calibrators: ['delay.K1','allcal_p.G0','allcal_p_jitter.G0','allcal_ap.G1_fluxscaled','bpcal.B0']
 - On targets: ['delay.K1','phscal_p_scan.G2','allcal_ap.G1_fluxscaled','bpcal.B0']


---
### 4.2.9 bandpass_1_sp

**bpcal_sp.B1**

Recalculate the BP table. Now the model columns contain the actual flux density of the calibrators, so the BP table will include spectral index information. I use `allcal_ap.G1_fluxscaled` table just because I want to have the right flux scale if applycal is selected for this step.

|Table parameter| Value                                                                    |
|-------------- | ------------------------------------------------------------------------ |
| name          | bpcal_sp.B1                                                              |
| field         | bpcal                                                                    |
| previous_cal  | ['delay.K1','allcal_p.G0','allcal_p_jitter.G0','allcal_ap.G1_fluxscaled']|
| solint        | inf                                                                      |
| spw           |                                                                          |
| uvrange       |                                                                          |
| combine       | field,scan                                                               |
| solnorm       | False                                                                    |
| table         | ./calib/inbase_bpcal_sp.B1                                               |
| gainfield     | bpcal                                                                    |
| interp        | nearest,linear                                                           |
| spwmap        | []                                                                       |

**Applycal:**

 - On calibrators: ['delay.K1','allcal_p.G0','allcal_p_jitter.G0','allcal_ap.G1_fluxscaled','bpcal_sp.B1']
 - On targets: ['delay.K1','phscal_p_scan.G2','allcal_ap.G1_fluxscaled','bpcal_sp.B1']


---
### 4.2.10 gain_1_amp_sp

**allcal_ap.G3**


|Table parameter| Value                                                        |
|---------------| ------------------------------------------------------------ |
| name          | allcal_ap.G3                                                 |
| field         | calsources                                                   |
| previous_cal  | ['delay.K1','allcal_p.G0','allcal_p_jitter.G0','bpcal_sp.B1']|
| solint        | 32s                                                          |
| gaintype      | G                                                            |
| calmode       | ap                                                           |
| spw           | *:innerchan                                                  |
| combine       |                                                              |
| table         | ./calib/inbase_allcal_ap.G3                                  |
| gainfield     | calsources                                                   |
| interp        | linear                                                       |
| spwmap        | []                                                           |

**allcal_ap_scan.G3**


|Table parameter| Value                                                        |
|---------------| ------------------------------------------------------------ |
| name          | allcal_ap_scan.G3                                            |
| field         | phscals                                                      |
| previous_cal  | ['delay.K1','allcal_p.G0','allcal_p_jitter.G0','bpcal_sp.B1']|
| solint        | inf                                                          |
| gaintype      | G                                                            |
| calmode       | ap                                                           |
| spw           | *:innerchan                                                  |
| combine       |                                                              |
| table         | ./calib/inbase_allcal_ap_scan.G3                             |
| gainfield     | phscals                                                      |
| interp        | linear                                                       |
| spwmap        | []                                                           |



**Applycal:**

 - On calibrators: ['delay.K1','bpcal_sp.B1','allcal_p.G0','allcal_p_jitter.G0','allcal_ap.G3']
 - On targets: ['delay.K1','bpcal_sp.B1','phscal_p_scan.G2','allcal_ap_scan.G3']

---
### 4.2.11 applycal_all
Two runs are executed. One to correct calibrators, using their own solutions when relevant. A second correction is executed for each target, using the solutions from the corresponding phase reference calibrator.

 - On calibrators: ['delay.K1','bpcal_sp.B1','allcal_p.G0','allcal_p_jitter.G0','allcal_ap.G3']
 - On targets: ['delay.K1','bpcal_sp.B1','phscal_p_scan.G2','allcal_ap_scan.G3']

---
### 4.2.12 flag_4_rflag

After calibration, we can run `flagdata` in rflag mode. This mode requires the data to be already calibrated.

CASA `flagdata`

|Parameter      | Value    |
|-------------- | -------- |
| mode          | rflag    |
| correlation   | ABS_ALL  |
| ntime         | 90min    |
| combinescans  | True     |
| datacolumn    | corrected|
| timedevscale  | 5        |
| freqdevscale  | 5        |
| action        | apply    |


---
### 4.2.13 plot_corrected

Produce plots of amp/phase vs time/freq for each baseline	plotms using corrected data column.

Inputs parameters needed:
```
inbase                  [str]
targets                 [str]
phscals                 [str]
fluxcal                 [str]
bpcal                   [str]
ptcal                   [str]
```

Output:
```
./plots/plots_corrected [directory with multiple plots in png format]
```

Produces plots in png format for the visibilities in `inbase_avg.ms`. It iterates through all sources specified by the user in the inputs file, and produces plots for each baseline. All plots are stored in the output directory. The plots have some averaging:

 - Amp/Phase vs Time: all channels averaged in each spw, colorized by spw.
 - Amp/Phase vs Frequency: averaged every 5 min, colorized by correlation.


---
### 4.2.14 weblog
Update the weblog will all available information and plots.

As the initial `summary_weblog` but the weblog will include all information related to the dataset, the observation, the calibration and the visibilities. There are four different web pages:

 - Home. Basic dataset information (name, date, antennas, frequency, averaging, etc.)
 - Observation summary. Access to listobs, antennas, elevation plot, uvcov plots per each source.
 - Calibration. List of tables produced and plots showing the solutions.
 - Plots. Amp/phase vs time/freq per source and per baseline. Corrected and uncorrected visibilities, and Amp/phase vs uvdist.

---

# 5. Support functions and variables

### get_msinfo [function] and msinfo [dict]
An internal task that retrieves information from the inputs file and also directly from the MS.

This is an inner function used by the pipeline, but the output can be useful. Produces the dictionary `msinfo` that contains: msfile, msfilename, project, run, sources, mssources, antennas, band, baselines, num_spw, t_ini, t_end, freq_ini, freq_end, chan_Res, nchan, innerchan, polarizations. The dictionary is saved with pickle in file `inbase.ms.msinfo.pkl` and `inbase_avg.ms.msinfo.pkl` if average data is produced. `msinfo` can be read with pickle, for example:

> `msinfo = pickle.load(open('inbase_avg.ms.msinfo.pkl', 'rb'))`


Example:

```
msfilename          : EGJ_3C264_L1_avg
nchan               : 128
run                 : EGJ_3C264_L1
t_ini               : 2015-02-19 18:07:14
innerchan           : 13~115
msfile              : ./EGJ_3C264_L1_avg.ms
polarizations       : R, L
chan_res            : 0.0005
num_spw             : 8
project             : Ex-gal J
freq_ini            : 1.2546495
band                : L
applycal_all        : True
antennas            : ['Cm' 'Da' 'De' 'Kn' 'Lo' 'Mk2' 'Pi']
freq_end            : 1.7661495
baselines           : ['Lo&Kn', 'Lo&De', 'Lo&Pi', 'Lo&Da', 'Lo&Cm', 'Mk2&Kn', 'Mk2&De', 'Mk2&Pi', 'Mk2&Da', 'Mk2&Cm', 'Kn&De', 'Kn&Pi', 'Kn&Da', 'Kn&Cm', 'De&Pi', 'De&Da', 'De&Cm', 'Pi&Da', 'Pi&Cm', 'Da&Cm']
t_end               : 2015-02-20 10:29:59
sources
   calsources       : 0319+415,1331+305,1407+284,1143+1834
   phscals          : 1143+1834
   no_fluxcal       : 0319+415,1407+284,1143+1834,1145+1936
   ptcal            : 0319+415
   cals_no_fluxcal  : 0319+415,1407+284,1143+1834
   bpcal            : 1407+284
   targets_phscals  : 1145+1936,1143+1834
   mssources        : 1145+1936,1407+284,1143+1834,1331+305,0319+415
   fluxcal          : 1331+305
   maincal          : 0319+415,1407+284,1331+305
   allsources       : 0319+415,1407+284,1143+1834,1145+1936,1331+305
   targets          : 1145+1936
   source_intent
      1145+1936     : targets
      1331+305      : fluxcal
      0319+415      : 
      1407+284      : bpcal
      1143+1834     : phscals

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





