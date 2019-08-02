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

1. [Optionally] Modify `default_params.json` or add manual flags to `inputfg_avg.flags` with your desired values
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
```

List of available steps to run/skip:

- pre_processing steps:
  - run_importfits  
  - flag_aoflagger  
  - flag_apriori        
  - flag_manual         
  - average             
  - plot_data           
  - save_flags          

- calibration steps:
  - restore_flags   
  - flag_manual_avg 
  - init_models         
  - bandpass            
  - initial_gaincal 
  - fluxscale           
  - bandpass_final  
  - gaincal_final   
  - applycal_all        
  - flag_target         
  - plot_corrected  
  - first_images    
  - split_fields  
  
**Examples of step selection**

You need to specify which steps of the pipeline to run. Some example on how to choose steps:

1. Run all the calibration steps (ideal for observatory-processed data for which you want to tweak the calibration parameters). Includes all calibrations steps (see list above):

`casa -c eMERLIN_CASA_pipeline/eMERLIN_CASA_pipeline.py -r calibration`

2. Run all pipeline steps:

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

Function `run_pipeline` parameters and defaults are: `run_pipeline(inputs_file='./inputs.ini', run_steps=[''], skip_steps=[''])`. Variables run_steps and skip_steps are python lists of steps as explained above.


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

**How do a fine tune the parameters used to run the tasks?**

All the parameters are included in the file `default_params.json` that you can edit manually. You can find the template inside the eMERLIN_CASA_pipeline folder. Steps to follow are:

- Copy the file `default_param.json` from the eMERLIN_CASA_pipeline directory to the working directory if it does not exist already.

```
cp eMERLIN_CASA_pipeline/default_params.json .
```

- Edit the file in the working directory and modify the appropriate parameters.
- Rerun the required steps of the pipeline. For example if you have just modified parameters for tasks in the calibration block, you can run:

```
casa -c eMERLIN_CASA_pipeline/eMERLIN_CASA_pipeline.py -r calibration
```

If you modify parameters affecting steps in the pre-processing block, you will need to rerun all the steps after the modified one.


You can have the file in your working directory where the pipeline is executed. If not there, the file inside the eMERLIN_CASA_pipeline directory will be used.

**I want to flag already processed data, how do I add manual flags?**
- Open/create a file called `inputfg_avg.flags` in the working directory.
- Write all your flag commands, one per line.
- You can find an example with the correct syntax in the [pipeline documentation](https://github.com/e-merlin/eMERLIN_CASA_pipeline/blob/master/documentation/docs.md#422-flag_manual_avg)
- Now you can execute the pipeline again selecting `-r calibration` (this will run the whole calibration part).
- Manual flags will be applied on the step `flag_manual_avg` of the calibration block.
- Verify that the flag commands file was correctly loaded, check the logs!

The pipeline accepts two optional manual flagging files:

- `inputfg.flags` will be applied to the unaveraged dataset in the step `flag_manual` (in the pre-processing block). Usually this file contains observatory flags and it is automatically created in the flag_apriori stage by the observatory staff.
- `inputfg_avg.flags` will be applied to the averaged dataset in the step `flag_manual_avg` (in the calibration block).

**How do I fill the source names in inputs.txt if I don't know which fields were observed?**

By default you should have all the information from the observatory. But if you only have the FITS-IDI and don't know the source names, you can run the first pipeline step alone `casa -c eMERLIN_CASA_pipeline/eMERLIN_CASA_pipeline.py -r run_importfits`. When the execution is finished, open the weblog and go to the tab `Observation Summary` where you will find the fields included in the MS and the listobs file with all the scans.

As a general rule, an observation will have 1331+3030 (3C286) as flux scale calibrator, 1407+2827 (OQ208) as bandpass calibrator and 0319+4130 (3C84) as bright ptcal calibrator. To distinguish between target and phasecal, you should look for alternating scans, and the target is usually the one with longer scans.
