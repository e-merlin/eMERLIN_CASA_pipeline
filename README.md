
## TL;DR ##
This is the CASA pipeline for e-MERLIN data, it is designed to be fully parallelised. It is, however, very early stages and will progress rapidly. It currently does work but does not do calibration yet. Feel free to change/modify.

## Dependencies ##
- aoflagger v2.7+ (v2.9+ recommended) (see https://sourceforge.net/projects/aoflagger/)
- CASA v4.6+ (see https://casa.nrao.edu/), tested on CASA v4.6.0
- python2.7

## Download ##
If you have git installed, you can get the pipeline using:  
`git clone https://github.com/e-merlin/CASA_eMERLIN_pipeline.git`

If you don't have git, you can download and unzip the files from [here](https://github.com/e-merlin/CASA_eMERLIN_pipeline/archive/master.zip).

## Usage ##
To run the pipeline simply do:
casa -c /path/to/pipeline/eMERLIN_CASA_pipeline.py -i <input file>`

To run the parallelized version using MPI in CASA you can use:
mpicasa -n <num_cores> casa -c /path/to/pipeline/eMERLIN_CASA_pipeline.py -i <input file>

#The pipeline uses a GUI interface as an input and now has a headless version. To run, do the following:
#   * Using GUI:  
#   `mpicasa -n <num_cores casa> -c /path/to/pipeline/eMERLIN_CASA_pipeline.py -g`
#   * With an input file (find sameple in pipeline. Define data location and steps to follow)  
#    `mpicasa -n <num_cores> casa -c /path/to/pipeline/eMERLIN_CASA_pipeline.py -i <input file>`
#  
#N.b. replace num_cores with the number of cores you would like to use.  
#To run not in parallel remove `mpicasa -n num_cores`

#IMPORTANT: The pipeline will look in your data directory for a ms or mms with the prefix defined in inputs and copy it to your current directory if it cannot find it in the current directory.


