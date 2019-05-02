This is the CASA pipeline for e-MERLIN data, it is designed to be fully parallelised. It is, however, very early stages and will progress rapidly. It currently does work but does not do calibration yet. Feel free to change/modify.

## Dependencies ##
- CASA v5.4+ (see https://casa.nrao.edu/)
- aoflagger v2.9+ (see https://sourceforge.net/projects/aoflagger/)

## Download ##
If you have git installed, you can get the pipeline using:  
`git clone https://github.com/e-merlin/eMERLIN_CASA_pipeline.git`

If you don't have git, you can download and unzip the files from [here](https://github.com/e-merlin/eMERLIN_CASA_pipeline/archive/master.zip).

To install other dependencies e.g. aoflagger/wsclean check out either A. Offringa's websites:
- aoflagger: https://sourceforge.net/projects/aoflagger/
- wsclean: https://sourceforge.net/projects/wsclean/

Or (recommended) use the handy anaconda scripts to instantly install dependcies within the conda environment. To do this follow the instructions in this repo.: https://github.com/jradcliffe5/radio_conda_recipes

To install other dependencies e.g. aoflagger/wsclean check out either A. Offringa's websites:
- aoflagger: https://sourceforge.net/projects/aoflagger/
- wsclean: https://sourceforge.net/projects/wsclean/

Or (**recommended!**) use the handy anaconda scripts to instantly install dependcies within the conda environment. To do this follow the instructions in this repo.: https://github.com/jradcliffe5/radio_conda_recipes

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



