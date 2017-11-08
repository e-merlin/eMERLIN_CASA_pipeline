<!---
Convert to pdf with http://www.markdowntopdf.com/
-->

# Inputs

There are two types of inputs. The ones at the begining expect a string that will depend on the project, sources, antennas or files needed. The ones in the last two blocks expect an integer that can be 0, 1 or 2 and are used to select which steps to execute. For `[str]` inputs, the use of single quotes is not required. Most inputs accept a list of values in the format of a comma-separated string. Examples: `target = '1111+2222'` or `target = 1111+2222`. For multiple inputs: `phscals = '1111+2222,3333+4444,5555+6666'`.


```
fits_path  [str]
```
Path to the location of the fits files. Can be an absolute or relative path. A single directory must be specified. All fits and FITS files from that directory are considered.


- `inbase     [str]`: Project name. It will be used to give name to MS, plots and tables.


- `targets    [str]`: Names of sources as they appear in the MS to be used as targets. Can be a comma-separated string. If more than one target is selected, the corresponding phase calibrator has to keep the same order.


- `phscals    [str]`: Names of sources as they appear in the MS to be used as calibrators. Can be a comma-separated string. If more than one phase calibrator is selected, the order needs to match the targets. If the same phase calibrator is used for different targets, just repeat the name as many times as needed.


- `fluxcal    [str]`: Names of source as it appears in the MS to be used as flux calibrator. Only one source accepted. 1331+305 expected. If a different source is selected the `fluxscale` step will not work properly.


- `bpcal      [str]`: Names of sources as they appear in the MS to be used as calibrators. Can be a comma-separated string.

- `ptcal      [str]`: Point-like calibrator. Names of sources as they appear in the MS to be included in calibration steps, but not used to calibrate other sources by now. For the moment, consider it for check sources.


- `refant     [str]`: Antenna name to be used as reference antenna. Can accept a comma-separated string, but CASA can only manage a single value for now (there are plans to allow `gaincal` to use a prioritized list. If empty, the pipeline will try to search for the most suitable antennas.


- `Lo_dropout_scans    [str]`: Comma-separated list of scans in which the Lo telescope was not observing the phase calibrator(s). This has two effects: those scans will be flagged in the `flag_1_apriori` step, and also each calibration table created will be edited to remove solutions for antenna Lo for the phase calibrator scans listed in `Lo_dropout_scans`. Example: `Lo_dropout_scans = '4,8,12,16,20,24,28'`

- `manual_flags_a      [str]`: Path to an external file that contains a list of flag commands, one per line. The flags will be applied to the unaveraged data set `inbase.ms`. The file will be read by CASA task flagdata using `mode='list'` and `inpfile` will be set to `manual_flags_a`. Note from the [CASA documentation](https://casa.nrao.edu/casadocs/casa-5.1.1/global-task-list/task_flagdata/about) on the format of the file: There should be no whitespace between KEY=VALUE since the parser first breaks command lines on whitespace, then on "=". Use only one whitespace to separate the parameters (no commas).

- `manual_flags_b      [str]`: Same as `manual_flags_a` but will be applied to the averaged data set `inbase_avg.ms`.





---
# Procedures

## 1. Pre-processing

### run_importfitsidi
Merge fits-IDI files in `fits_path` to form a MS named `inbase.ms`

Inputs parameters needed:
```
fits_path               [str]
inbase                  [str]
```

Produces:
```
inbase.ms               [MS]
inbase.ms.listobs       [txt]
inbase.ms.sp0           [MS] (only for mixed mode observations)
```

First `inbase.ms` is removed if present (`inbase` is the parameter in inputs file). All files ending in .fits or .FITS in the folder `fits_path` are considered for importing and will be merged. The logger shows which files have been found. Several CASA tasks are executed: `importfitsidi` is executed with `constobsid=True, scanreindexgap_s=15.0`. Then `uvfix` with option `reuse=False` is run to correct UVW values and the output MS replaces `inbase.ms`. Then `flagdata` with `mode='manual',autocorr=True` to remove autocorrelations. Finally, `listobs` is run and a new listobs file is created in `inbase.ms.listobs`.

At this point the pipeline checks if the observations were observed in mixed mode. If so, the spectal spws will be splitted and `inbase.ms` will contain only the continuum data (broadband, low spectral resolution). The high resolution spw will be splitted in new MS ending with sp0, sp1, etc.


---
### summary_weblog
Retrieve basic information from the original MS, produce some plots and compile eveything in the weblog.

Inputs parameters needed:
```
inbase                  [str]
```

Produces:
```
inbase.ms.msinfo.pkl    [pickle dictionary]
```
Runs `get_msinfo` (see section 3, below) to get information from the inputs file and directly from the MS. It prepares plots for elevation vs time and uvcov for each individual sources. It produces a weblog with the available information. If no `refant` is selected, or the antenna selected is not in the MS, a procedure will try to guess a good refant to use. Please update the inputs file with a suitable refant.

---
### hanning
Just runs CASA hanning smoothing.

Inputs parameters needed:
```
inbase                  [str]
```

Produces:
```
inbase.ms               [MS]
```

Runs mstransform with mode='hanning' on DATA column. It produces a new MS, which substitutes the `inbase.ms`. It works with MS or MMS. This step is optional, but recommended for L band observations. It is probably not needed for C band datasets.



---
### ms2mms
Optional step to transform original MS to a MMS.

Inputs parameters needed:
```
inbase                  [str]
```

Produces:
```
inbase.mms               [MMS]
```

This allow CASA to run tasks on the MultiMeasurementSet in parallel transparently. Takes time to convert but later makes calibration and imaging faster. To be tested in detail. Uses CASA task `partition` with default parameters. **Important**: if ms2mms is executed, it will delete `inbase.ms` and all the following steps will be run on the new `inbase.mms`.


---
### flag_0_aoflagger
Uses aoflagger to autoflag data using predefined strategies

Inputs parameters needed:
```
inbase                  [str]
```

It runs [aoflagger](https://sourceforge.net/p/aoflagger/wiki/Home/) on all the fields in the MS. aoflagger v2.9+ required. The pipeline might work with previous versions in some situations, but it is not recommended and probably will fail The strategies are selected by field following this criteria: (1) First, check if user has a new strategy for this field in the local folder ./aoflagger_strategies/user/, the strategy needs to match: `<field>.rfis`. (2) If not, check if user has produced a new strategy for this field in the pipeline folder. So will search for `pipeline_path+'aoflagger_strategies/default/<field>.rfis`, (3) If nothing is found, just use the default strategy. For example for field 1234+5678 the prioritization for searching strategies is:

 1. ./aoflagger_strategies/user/1234+5678.rfis
 2. /path/to/pipeline/aoflagger_strategies/default/1234+5678.rfis
 3. /path/to/pipeline/aoflagger_strategies/default/default_faint.rfis



---
### flag_1_apriori
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
 - Quack 20 seconds for all targets and phasecals.
 - Scans in `Lo_dropout_scans` for antenna = 'Lo' and fields in `phscals`.


---
### flag_2a_manual
Applies flags from external file with a list of flag commands

```
inbase                  [str]
manual_flags_a          [str, path to file]
```

Simply runs `flagdata(vis='inbase.ms', mode='list', inpfile=manual_flags_a)`. For more information on this flagdata mode see [flagdata](https://casa.nrao.edu/casadocs/casa-5.1.1/global-task-list/task_flagdata/about). For details on the format of the file see [examples](https://casa.nrao.edu/casadocs/casa-5.1.1/global-task-list/task_flagdata/examples). It needs one flag command per line, with parameters separated by one space (don't use commas!). For example:

```
mode='manual' field='1331+305' antenna='' timerange='10:00:00~10:11:30'
mode='manual' field='' antenna='' timerange='' spw='0:0~30'
mode='manual' field='' antenna='Mk2' timerange='09:05:00~16:27:00'
mode='manual' field='1258-2219' antenna='' timerange='12:57:01~12:59:59'
mode='quack' field='1258-2219,1309-2322' quackinterval=24.
```
---
### average_1
Split dataset and average to reduce data volume


---
### plot_data
Produce plots of amp/phase vs time/freq for each baseline	plotms




## 2. Calibration

### flag_2b_manual

### init_models

### bandpass_0

### flag_3_tfcropBP

### delay

### gain_0_p_ap

### fluxscale

### bandpass_1_sp

### gain_1_amp_sp

### applycal_all

### flag_4_rflag

### plot_corrected

### weblog


## 3. Support functions


### get_msinfo
Internal task that retrieves information from the inputs file and also directly from the MS.

This is an inner funtion used by the pipeline, but the output can be useful. Produces the dictionary `msinfo` that contains: msfile, msfilename, project, run, sources, mssources, antennas, band, baselines, num_spw, t_ini, t_end, freq_ini, freq_end, chan_Res, nchan, innerchan, polarizations. The dictionary is saved with pickle in file `inbase.ms.msinfo.pkl` and `inbase_avg.ms.msinfo.pkl` if average data is produced. `msinfo` can be read with pickle or with:

> `msinfo = pickle.load('inbase.ms.msinfo.pkl')`
> 
> `em.prt_hist(msinfo)` # To see the contents

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


---
# Quick summary

### Pre-processing
|Procedure        |Summary                                                     | CASA/external tasks                       |Inputs                   |Outputs        |Notes                                                                                          |
|-----------------|------------------------------------------------------------|-------------------------------------------|-------------------------|---------------|-----------------------------------------------------------------------------------------------|
|run_importfitsidi| Concatenate all fits-IDI files in a folder to a MS         | importfitsidi, fixvis, flagdata (autocorr)|fits_path                |inbase.ms      | Also produces inbase.ms.listobs                                                               |
|summary_weblog   | Produces weblog with basic project information and plots   | get_data, plotms                          |                         |               | Info (sources, spw, channels, time, etc), elevation vs time plot and uvcov plots.             |
| hanning         | Run hanning smoothing                                      | mstransform                               |                         |               | Hanning smoothing helps with RFI removal. Recommended for L band data.                        |                                                                                                                                                                              | ms2mms          | Convert to mms format to run tasks in parallel             | mstransform                               |                         |inbase.mms     | If this step is executed, inbase.ms will be removed and all other steps will run on inbase.mms|
| flag0_aoflagger | Autoflag data using predefined strategies                  | aoflagger                                 |                         |               | Accepts user strategies per field. Can take very long if low memory available                 |
| flag1_apriori   | Standard flags: Lo-Mk2, edge channels, slewing, Lo dropouts| flagdata                                  |sources, Lo_dropout_scans|               |                                                                                               |
| flag_2a_manual  | Applies flags from external file with a list of flag commands | flagdata                               |manual_flags_a           |               | Applies flag commands prepared by the user in external file to **unaveraged** data            |
| average_1       | Split dataset and average to reduce data volume            | split                                     |sources                  |inbase_avg.ms  |                                                                                               |
| plot_data       | Produce plots of amp/phase vs time/freq for each baseline  | plotms                                    |                         |plots (several)| Plots DATA column only.                                                                       |



### Calibration

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





