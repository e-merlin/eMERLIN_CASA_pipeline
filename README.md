1. [ Description ](#description)
1. [ Dependencies ](#dependencies)
1. [ Download ](#download)
1. [ Quick start ](#quickstart)
1. [ Usage ](#usage)
1. [ Additional information ](#information)
1. [ FAQ ](#faq)



<a name="description"></a>
## Description ##

The e-MERLIN CASA Pipeline (eMCP) is a python pipeline working on top of [CASA](https://casa.nrao.edu/) to process and calibrate interferometric data from the [e-MERLIN](http://www.e-merlin.ac.uk/) array. Access to data information, statistics and assessment plots on calibration tables and visibilities can be accessed by the pipeline weblog, which is updated in real time as the pipeline job progresses. The output is calibrated data and preliminary lookup images of the relevant fields. It can calibrate mixed mode data that includes narrow-band high spectral resolution spectral windows for spectral lines, and also special observing modes as pseudo-wideband observations. Currently no polarization calibration is performed.

<a name="dependencies"></a>
## Dependencies ##
- CASA v5.5+ (see https://casa.nrao.edu/)
- aoflagger v2.9+ (see https://sourceforge.net/projects/aoflagger/), only needed to calibrate L-band data.

<a name="download"></a>
## Download ##
If you have git installed, you can get the pipeline using:  
`git clone https://github.com/e-merlin/eMERLIN_CASA_pipeline.git`

If you don't have git, you can download and unzip the files from [here](https://github.com/e-merlin/eMERLIN_CASA_pipeline/archive/master.zip).

To install aoflagger check out either A. Offringa's websites:
- aoflagger: https://sourceforge.net/projects/aoflagger/
- wsclean: https://sourceforge.net/projects/wsclean/  (not required for running the pipeline)

or (recommended) use the handy anaconda scripts to instantly install dependcies within the conda environment. To do this follow the instructions in this repo.: https://github.com/jradcliffe5/radio_conda_recipes

<a name="quickstart"></a>
## Quick start ##

If you have received calibrated data from the observatory and you want to refine the calibration, you can:

1. [Optionally] Modify `default_params.json` or add manual flags to `manual_avg.flags` with your desired values.
2. Run:

`casa -c eMERLIN_CASA_pipeline/eMERLIN_CASA_pipeline.py -r calibration`

<a name="usage"></a>
## Usage ##
Normal pipeline execution. When you have in your working directory the file `inputs.ini` and you have extracted the pipeline:

`casa -c /path/to/pipeline/eMERLIN_CASA_pipeline.py`

To run the parallelized version using MPI in CASA you can use:  

`mpicasa -n <num_cores> casa -c /path/to/pipeline/eMERLIN_CASA_pipeline.py`

**Optional arguments**

Names in capital need to be set by the user:

```
  -h, --help                     show this help message and exit
  
  
  -i INPUTS_FILE
  --inputs INPUTS_FILE
                                 Inputs file to use. Default is inputs.ini
                                 
                                 
  -r          RUN_STEPS [RUN_STEPS ...]
  --run-steps RUN_STEPS [RUN_STEPS ...]
                                 Whitespace separated list of steps to run. Apart from
                                 individual steps, it also accepts "all",
                                 "pre_processing" and "calibration"
                                 
                                 
  -s           SKIP_STEPS [SKIP_STEPS ...]
  --skip-steps SKIP_STEPS [SKIP_STEPS ...]
                                 Whispace separated list of steps to skip
                                 
                                 
  -l
  --list-steps                   Show list of available steps and exit

```

You can get the list of available steps with:

`casa -c eMERLIN_CASA_pipeline/eMERLIN_CASA_pipeline.py -l`

```
pre_processing
    run_importfits
    flag_aoflagger
    flag_apriori
    flag_manual
    average
    plot_data
    save_flags
    
calibration
    restore_flags
    flag_manual_avg
    init_models
    bandpass
    initial_gaincal
    fluxscale
    bandpass_final
    gaincal_final
    applycal_all
    flag_target
    plot_corrected
    first_images
    split_fields

```

Selection options are any combination of: a list of any individual step names, `pre_processing`, `calibration` or `all`
  
**Examples of step selection**

You need to specify which steps of the pipeline to run. Some example on how to choose steps:

1. Run all the calibration steps (ideal for observatory-processed data for which you want to tweak the calibration parameters). Includes all calibrations steps (see list above):

`casa -c eMERLIN_CASA_pipeline/eMERLIN_CASA_pipeline.py -r calibration`

2. Run all pipeline steps (you will need the raw FITS-IDI files for the initial step):

`casa -c eMERLIN_CASA_pipeline/eMERLIN_CASA_pipeline.py -r all`

3. Run only the pre-processing steps (usually executed by the observatory. Otherwise you need the raw FITS-IDI files):

`casa -c eMERLIN_CASA_pipeline/eMERLIN_CASA_pipeline.py -r pre_processing`

4. Any combination of the steps above, for example:

`casa -c eMERLIN_CASA_pipeline/eMERLIN_CASA_pipeline.py -r plot_corrected first_images split_fields`

5. Run all calibration steps except plot_corrected:

`casa -c eMERLIN_CASA_pipeline/eMERLIN_CASA_pipeline.py -r calibration -s plot_corrected`

**Running the pipeline interactively from CASA**

To execute the pipeline from a running CASA instance you need to write in the CASA shell:

~~~~
run_in_casa = True
pipeline_path = '/path/to/pipeline_path/'   # You need to define this variable explicitly
execfile(pipeline_path + 'eMERLIN_CASA_pipeline.py')
eMCP = run_pipeline(run_steps=['calibration'])
~~~~

Function `run_pipeline` parameters and defaults are: `run_pipeline(inputs_file='./inputs.ini', run_steps=[], skip_steps=[])`. Variables run_steps and skip_steps are python lists of steps as explained above.


<a name="information"></a>
## Additional information ##

- [Documentation [online]](documentation/docs.md)
- [Wiki pages](https://github.com/e-merlin/eMERLIN_CASA_pipeline/wiki)


<a name="faq"></a>
## FAQ ##

**How do I open the weblog?**

The weblog consist of a series of html files. From the working directory you can open the file `./weblog/index.html` with your preferred web browser.

**How do I know what has been executed?**

You can visit the tab `Pipeline info` in the weblog, where you will find which steps were executed. You will also find a link to the Pipeline log, the CASA log and two files with all the parameters used during the data processing.

**I want to re-run the pipeline to improve the calibration, what do I change?**

There are two main blocks: pre-processing and calibration. Most probably you will only need to repeat the calibration part. Recommended course of action:

- Identify changes you want to include in the data reduction, like changing calibration parameters or adding manual flags.
- Add or edit file `manual_avg.flags` with your flag commands (follow the CASA syntax).
- Edit the file `inputs.ini` if you need to change the sources used or they intend.
- Edit the file `default_params.json` changing any parameter the pipeline is using, if needed.
- Run the calibration block of the pipeline with the command:

`casa -c ./eMERLIN_CASA_pipeline/eMERLIN_CASA_pipeline.py -r calibration`

**Which flag files does the pipeline accept and what is the right syntax?**

There are four different flag files accepted by the pipeline:

| Flag file     | Used by step  | Notes |
| ------------- |:-------------:| -----|
| observatory.flags | flag_apriori | Created by the observatory with antenna slewing or other major faults. Please donot edit it yourself. |
| manual.flags | flag_manual | This is meant to flag the unaveraged data set during the pre-processing stage |
| manual_avg.flags | flag_manual_avg | This is meant to flag the averaged data set during the calibration stage |
| manual_narrow.flags | flag_manual_avg | Use this to add flag commands for narrow-band spectral line data set|

For the syntax needed for CASA follow [Basic Syntax Rules](https://casa.nrao.edu/casadocs/casa-5.5.0/global-task-list/task_flagdata/about) in the CASA docs flagdata (end of the section). The main rules are:

1. Use only ONE white space to separate the parameters (no commas). Each key should only appear once on a given command line/string.
2. There is an implicit mode for each command, with the default being 'manual' if not given.
3. Comment lines can start with '#' and will be ignored. The parser used in flagdata will check each parameter name and type and exit with an error if the parameter is not a valid flagdata parameter or of a wrong type.

Example for e-MERLIN:

```
mode='manual' field='1331+305' antenna='' timerange='10:00:00~10:11:30'
mode='manual' field='' antenna='' timerange='' spw='0:0~30'
mode='manual' field='' antenna='Mk2' timerange='09:05:00~16:27:00'
mode='manual' field='1258-2219' antenna='' timerange='12:57:01~12:59:59'
mode='quack' field='1258-2219,1309-2322' quackinterval=24.
```

Example from the CASA docs:
```
scan='1~3' mode='manual'
# this line will be ignored
spw='9' mode='tfcrop' correlation='ABS_XX,YY' ntime=51.0
mode='extend' extendpols=True
scan='1~3,10~12' mode='quack' quackinterval=1.0
```

**How do I fill the source names in inputs.ini if I don't know which fields were observed?**

By default you should have all the information from the observatory. But if you only have the FITS-IDI and don't know the source names, you can run the first pipeline step alone `casa -c eMERLIN_CASA_pipeline/eMERLIN_CASA_pipeline.py -r run_importfits`. When the execution is finished, open the weblog and go to the tab `Observation Summary` where you will find the fields included in the MS and the listobs file with all the scans.

As a general rule, an observation will have 1331+3030 (3C286) as flux scale calibrator, 1407+2827 (OQ208) as bandpass calibrator and 0319+4130 (3C84) as bright ptcal calibrator. To distinguish between target and phasecal, you should look for alternating scans, and the target is usually the one with longer scans.
