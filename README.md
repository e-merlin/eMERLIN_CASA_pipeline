The e-MERLIN CASA Pipeline (eMCP) is a python pipeline working on top of [CASA](https://casa.nrao.edu/) to process and calibrate interferometric data from the [e-MERLIN](http://www.e-merlin.ac.uk/) array. Access to data information, statistics and assessment plots on calibration tables and visibilities can be accessed by the pipeline weblog, which is updated in real time as the pipeline job progresses. The output is calibrated data and preliminary lookup images of the relevant fields. It can calibrate mixed mode data that includes narrow-band high spectral resolution spectral windows for spectral lines, and also special observing modes as pseudo-wideband observations. Currently no polarization calibration is performed.


## Dependencies ##
- CASA v5.4+ (see https://casa.nrao.edu/)
- aoflagger v2.9+ (see https://sourceforge.net/projects/aoflagger/), only needed to calibrate L-band data.

## Download ##
If you have git installed, you can get the pipeline using:  
`git clone https://github.com/e-merlin/eMERLIN_CASA_pipeline.git`

If you don't have git, you can download and unzip the files from [here](https://github.com/e-merlin/eMERLIN_CASA_pipeline/archive/master.zip).

To install other dependencies e.g. aoflagger/wsclean check out either A. Offringa's websites:
- aoflagger: https://sourceforge.net/projects/aoflagger/
- wsclean: https://sourceforge.net/projects/wsclean/

Or (recommended) use the handy anaconda scripts to instantly install dependcies within the conda environment. To do this follow the instructions in this repo.: https://github.com/jradcliffe5/radio_conda_recipes

## Usage ##
The easiest case, when you have in your working directory the file inputs.txt and you have extracted the pipeline:

`casa -c eMERLIN_CASA_pipeline/eMERLIN_CASA_pipeline.py -i inputs.txt`

A more general case would be:

`casa -c /path/to/pipeline/eMERLIN_CASA_pipeline.py -i <input file>`

To run the parallelized version using MPI in CASA you can use:  

`mpicasa -n <num_cores> casa -c /path/to/pipeline/eMERLIN_CASA_pipeline.py -i <input file>`

To execute the pipeline from within CASA:
~~~~
> run_in_casa = True
> pipeline_path = '/path/to/pipeline_path/'   # You need to define this variable explicitly
> execfile(pipeline_path + 'eMERLIN_CASA_pipeline.py')
> inputs, msinfo = run_pipeline(inputs_path=<input file>)
~~~~

## Additional information ##

- [Documentation [online]](documentation/docs.md)
- [Wiki pages](https://github.com/e-merlin/eMERLIN_CASA_pipeline/wiki)



## FAQ ##

### How do I open the weblog?

The weblog consist of a series of html files. From the working directory you can open the file `./weblog/index.html` with your preferred web browser.

### How do I fill the source names in inputs.txt if I don't know which fields were observed?

By default you should have all the information from the observatory. But if you only have the FITS-IDI and don't know the source names, you can run the first pipeline step alone `run_importfits = 1` and set all other steps to `0`. When the execution is finished, open the weblog and go to the tab `Observation Summary` where you will find the fields included in the MS and the listobs file with all the scans.

As a general rule, an observation will have 1331+3030 (3C286) as flux scale calibrator, 1407+2827 (OQ208) as bandpass calibrator and 0319+4130 (3C84) as bright ptcal calibrator. To distinguish between target and phasecal, you should look for alternating scans, and the target is usually the one with longer scans.

### How do I know what has been executed?

You can visit the tab `Pipeline info` in the weblog, where you will find which steps were executed. You will also find a link to the Pipeline log, the CASA log and two files with all the parameters used during the data processing.

### How do a fine tune the parameters used to run the tasks?

All the parameters are included in the file `default_params.json` that you can edit manually. You can have the file in your working directory where the pipeline is executed. If not there, the file inside the eMERLIN_CASA_pipeline directory will be used.
