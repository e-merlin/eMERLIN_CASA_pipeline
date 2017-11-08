This is the CASA pipeline for e-MERLIN data, it is designed to be fully parallelised. It is, however, very early stages and will progress rapidly. It currently does work but does not do calibration yet. Feel free to change/modify.

## Dependencies ##
- CASA v5.0+ (see https://casa.nrao.edu/)
- aoflagger v2.9+ (see https://sourceforge.net/projects/aoflagger/). v2.7+ accepted but not recommended, will be deprecated.
- python2.7

## Download ##
If you have git installed, you can get the pipeline using:  
`git clone https://github.com/e-merlin/CASA_eMERLIN_pipeline.git`

If you don't have git, you can download and unzip the files from [here](https://github.com/e-merlin/CASA_eMERLIN_pipeline/archive/master.zip).

## Usage ##
To run the pipeline simply do:  
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

- [Documentation [markdown]](documentation/docs.md)
- [Documentation [html]](documentation/docs.html)
- [Wiki pages](https://github.com/e-merlin/CASA_eMERLIN_pipeline/wiki)


